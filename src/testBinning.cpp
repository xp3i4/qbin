// ==========================================================================
//                           Mapping SMRT reads 
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: cxpan <chenxu.pan@fu-berlin.de>
// ==========================================================================

#include <csignal>
#include "mapper.h"

using namespace seqan; 

template <typename TDna, typename TSpec>
Mapper<TDna, TSpec>::Mapper(Options & options):
    record(options),
    qIndex(genomes()),
    of(toCString(options.getOutputPath()))
{
        switch (options.sensitivity)
        {
            case 0: 
            {
                parm = parm0; //normal
                break;
            }
            case 1:
            {
                parm =  parm1; //fast
                break;
            }
            case 2:
            {
                parm = parm2; //sensitive
                break;
            }
        }
        _thread = options.thread;
        
        std::cerr << "[mapper thread] " << _thread << "\n";
        
}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex()
{
    std::cerr << ">[Creating index] \n";
    float ythredfrac = 0.8;
    createHIndex(genomes(), bin(), qIndex, ythredfrac, _thread);
    return 0;
}

template <typename TDna, typename TSpec>
int Mapper<TDna, TSpec>::createIndex2_MF()
{
    std::cerr << ">[Creating index] \n";
    float ythredfrac = 0.8;
    createHIndex2_MF(genomes(), bin(), qIndex, ythredfrac, _thread);
    return 0;
}


template <typename TDna, typename TSpec>
inline unsigned testbin(typename PMCore<TDna, TSpec>::Index & index,
                        typename PMRecord<TDna>::RecSeqs & reads,
                        StringSet<String<uint64_t> > & list,
//                        MapParm & mapParm,
                        unsigned binNo,
                        unsigned threads
                             )
{   
    std::cerr << "[degbu]::binNO "<< binNo << "\n";
    double time = sysTime();
    unsigned step = 1;
    typedef typename PMCore<TDna, TSpec>::Index TIndex;
    typedef typename TIndex::TShape PShape;
    unsigned ysthred = 0;
    //std::cerr << "[debug] " << threads << "\n";
#pragma omp parallel
{
    PShape shape;
    unsigned dt = 0;
    String<unsigned> score; //(binNo);
    resize (score, binNo, 0);
    StringSet<String<uint64_t> > tmpRslt;
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;

    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(tmpRslt, ChunkSize);
    unsigned c = 0;
    
#pragma omp for 
    for (unsigned j = 0; j < length(reads); j++)
    {
        hashInit(shape, begin(reads[j]));
        for (unsigned k = 0; k < length(reads[j]) - shape.span + 1; k++)
        {
            hashNext(shape, begin(reads[j]) + k);

            if (++dt == step)
            {
                uint64_t pos = getXDir(index, shape.XValue, shape.YValue);
                while (_DefaultHs.isBody(index.ysa[pos]))
                {
                    if (_DefaultHs.getHsBodyY(index.ysa[pos]) == shape.YValue)
                    {
                        score[_DefaultHs.getHsBodyS(index.ysa[pos])] += 1;
                    }
                    ++pos;
                }
                dt = 0;
            }
            
        }
        for (unsigned k = 0; k < length(score); k++)
        {
            if (score[k] > ysthred)
            {
                appendValue(tmpRslt[c], k);
            }
            score[k] = 0;
        }
        c += 1;
    }
#pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
#pragma omp ordered
    {
        append(list, tmpRslt);
    }
}
std::cerr << ">mapping[s] " << sysTime() - time << "\n";
    return 0;
}


template <typename TDna, typename TSpec>
int map(Mapper<TDna, TSpec> & mapper)
{
    //printStatus();
    omp_set_num_threads(mapper.thread());
    //mapper.createIndex(); // true for parallel 
    mapper.createIndex(); 
    SeqFileIn rFile(toCString(mapper.readPath()));
    
    
    double time = sysTime();
    std::cerr <<  ">reading reads from " << mapper.readPath() << "\r";
    readRecords(mapper.readsId(), mapper.reads(), rFile);//, blockSize);
    std::cerr << ">end reading " <<sysTime() - time << "[s]" << std::endl;
    std::cerr << ">mapping " << length(mapper.reads()) << " reads to reference genomes"<< std::endl;
    testbin<TDna, TSpec>(mapper.index(), mapper.reads(), mapper.rslt(), length(mapper.bin()), mapper.thread());
    
    std::cerr << ">writing result to disk \n";
    for (unsigned k = 0; k < length(mapper.rslt()); k++)
    {
        mapper.of_stream() << "read_" << k << " ";
        for (unsigned j = 0; j < length(mapper.rslt()[k]); j++)
        {
            mapper.of_stream() << mapper.rslt()[k][j] << " ";
        }
        mapper.of_stream() << "\n";
    }

    return 0;
}

seqan::ArgumentParser::ParseResult
parseCommandLine(Options & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("pacMapper");
    // Set short description, version, and date.
    setShortDescription(parser, "Alignment of SMRT sequencing read");
    setVersion(parser, "1.0");
    setDate(parser, "May 2017");

    // Define usage line and long description.
    addUsageLine(parser,
                    "[\\fIOPTIONS\\fP] \"\\fIread.fa\\fP\" \"\\fIgnome.fa\\fP\"");
    addDescription(parser,
                    "Program for mapping raw SMRT sequencing reads to reference genome.");

    // Argument.
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "read"));
    setHelpText(parser, 0, "Reads file .fa, .fasta");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "genome", true));
    setHelpText(parser, 1, "Reference file .fa, .fasta");

    addSection(parser, "Mapping Options");
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "choose output file.",
            seqan::ArgParseArgument::STRING, "STR"));
    addOption(parser, seqan::ArgParseOption(
        "s", "sensitivity", "Sensitivity mode. -s 0 normal {DEFAULT} -s 1 fast  -s 2 sensitive",
            seqan::ArgParseArgument::INTEGER, "INT"));
    addOption(parser, seqan::ArgParseOption(
        "t", "thread", "Default -t 4",
            seqan::ArgParseArgument::INTEGER, "INT"));
        
    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser,
                "\\fBpacMapper\\fP \\fB-U\\fP \\fIchr1.fa reads.fa\\fP",
                "Print version of \"rd\"");
    addListItem(parser,
                "\\fBpacMapper\\fP \\fB-L\\fP \\fB-i\\fP \\fI3\\fP "
                "\\fIchr1.fa reads.fa\\fP",
                "Print \"\" with every third character "
                "converted to upper case.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    getOptionValue(options.oPath, parser, "output");
    getOptionValue(options.sensitivity, parser, "sensitivity");
    getOptionValue(options.thread, parser, "thread");

    seqan::getArgumentValue(options.rPath, parser, 0);
    options.gPath = seqan::getArgumentValues(parser, 1);
    //for (unsigned k = 0; k < length(options.gPath); k++)
    //    std::cout << "[debug]::g " << " " << options.gPath[k] << std::endl;

    return seqan::ArgumentParser::PARSE_OK;

}

int main(int argc, char const ** argv)
{
    double time = sysTime();

    std::cerr << "Encapsulated version: Binning by q-gram index" << std::endl;
    (void)argc;
    // Parse the command line.
    Options options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    Mapper<> mapper(options);
    //mapper.printParm();
    //std::cout << "[debug]::genomePath " << mapper.genomePath() << std::endl;
    map(mapper);
    std::cerr << "Time in sum[s] " << sysTime() - time << std::endl;

    return 0;
}
