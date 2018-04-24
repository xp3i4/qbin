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

#include "base.h"
//#include "pmpfinder.h"
#include "mapparm.h"

#ifndef SEQAN_HEADER_PACMAPPER_H
#define SEQAN_HEADER_PACMAPPER_H

template <typename TDna = typename MapperBase<>::DefaultAlphabet, 
    typename TSpec = typename MapperBase<>::DefaultShape>
class Mapper {
    typedef MapperBase<TDna, TSpec> Base;
    typedef typename Base::MRecord   Record;
    typedef typename Base::MParm     Parm;
    typedef typename Base::MIndex    Index;
    typedef typename Base::MAnchors    Anchors;
    typedef typename Base::MRes   Res;
    typedef typename Base::MSeq      Seq;
    typedef typename Base::MSeqs     Seqs;
    typedef typename Res::HitSet    HitSet;
    typedef typename Res::HitType   HitType; 
    typedef StringSet<String<uint64_t> > Rst;

    Record  record;
    Parm    parm;
    Res     res;
    Index   qIndex;
    std::ofstream of;
    unsigned _thread;
    Rst rst;

public:
    Mapper();
    Mapper(Options & options);
    Seqs & reads() {return record.seq1;}             
    Seqs & genomes() {return record.seq2;}             
    Parm & mapParm() {return parm;}
    Res & result() {return res;}
    Index & index() {return qIndex;}
    
    void printHits();
    void printBestHitsStart();
    void printResult();
    void printParm();
    int createIndex();
    int createIndex2_MF();//destruct genomes string during the creation to reduce memory footprint
    unsigned sens(){return parm.sensitivity;}
    unsigned & thread(){return _thread;}
    CharString & readPath(){return record.readPath;}
    CharString & genomePath(){return record.genomePath;}
    StringSet<CharString> & readsId(){return record.id1;}
    StringSet<CharString> & genomesId(){return record.id2;}
    String<uint64_t>  & bin(){return record.bin;}
    Rst & rslt(){return rst;}
    std::ofstream & of_stream() {return of;}
    
};



#endif
