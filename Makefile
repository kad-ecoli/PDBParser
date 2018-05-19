CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=qTMclust pdb2fasta reindex_pdb strip_sidechain split_chain BackboneTorsion SidechainTorsion NWalign pdb2rmsd pdb2ss pdb2sarst PeptideBuilder addCbeta
OLD_PROG=pdb2ThreeDblast SarstAlign ThreeDblastAlign ccealign getRotSeq
HEADER=PDBParser.hpp pstream.h
NW_HEADER=NWalign.hpp FilePathParser.hpp ROTSUMalign.hpp getRotSeq.hpp SidechainTorsion.hpp GeometryTools.hpp MathTools.hpp SarstAlign.hpp StructuralAlphabet.hpp BackboneTorsion.hpp ThreeDblastAlign.hpp SSalign.hpp 

current: ${PROG}

all: ${PROG} ${OLD_PROG}

pdb2fasta: pdb2fasta.cpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

split_chain: split_chain.cpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

reindex_pdb: reindex_pdb.cpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

strip_sidechain: strip_sidechain.cpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

BackboneTorsion: BackboneTorsion.cpp BackboneTorsion.hpp GeometryTools.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

SidechainTorsion: SidechainTorsion.cpp SidechainTorsion.hpp GeometryTools.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

getRotSeq: getRotSeq.cpp getRotSeq.hpp SidechainTorsion.hpp GeometryTools.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2sarst: pdb2sarst.cpp StructuralAlphabet.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2ThreeDblast: pdb2ThreeDblast.cpp StructuralAlphabet.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2ss: pdb2ss.cpp StructuralAlphabet.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp ${NW_HEADER} ${HEADER} pdb2rmsd.hpp TMalign.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS} -ITMalign

SarstAlign: SarstAlign.cpp SarstAlign.hpp StructuralAlphabet.hpp NWalign.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

ThreeDblastAlign: ThreeDblastAlign.cpp ThreeDblastAlign.hpp StructuralAlphabet.hpp NWalign.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2rmsd: pdb2rmsd.cpp pdb2rmsd.hpp NWalign.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

qTMclust: qTMclust.cpp qTMclust.hpp ${NW_HEADER} ${HEADER} pdb2rmsd.hpp TMalign.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS} -ITMalign

ccealign: ccealign.cpp ccealign.hpp ce2tm.hpp ${NW_HEADER} ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS} -Itnt

PeptideBuilder: PeptideBuilder.cpp PeptideBuilder.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

addCbeta: addCbeta.cpp PeptideBuilder.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm ${PROG}
	rm ${OLD_PROG}
