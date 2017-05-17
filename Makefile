CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=pdb2fasta reindex_pdb strip_sidechain split_chain BackboneTorsion SidechainTorsion NWalign pdb2rmsd
OLD_PROG=pdb2sarst pdb2ThreeDblast SarstAlign ThreeDblastAlign
HEADER=PDBParser.hpp pstream.h

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

pdb2sarst: pdb2sarst.cpp StructuralAlphabet.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2ThreeDblast: pdb2ThreeDblast.cpp StructuralAlphabet.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp NWalign.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

SarstAlign: SarstAlign.cpp SarstAlign.hpp StructuralAlphabet.hpp NWalign.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

ThreeDblastAlign: ThreeDblastAlign.cpp ThreeDblastAlign.hpp StructuralAlphabet.hpp NWalign.hpp BackboneTorsion.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdb2rmsd: pdb2rmsd.cpp pdb2rmsd.hpp NWalign.hpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm ${PROG}
	rm ${OLD_PROG}
