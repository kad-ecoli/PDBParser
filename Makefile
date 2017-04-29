CC=g++
CFLAGS=-O3
LDFLAGS=-static

PROG=pdb2fasta reindex_pdb strip_sidechain split_chain BackboneTorsion
HEADER=PDBParser.hpp pstream.h GeometryTools.hpp

all: ${PROG}

pdb2fasta: pdb2fasta.cpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

split_chain: split_chain.cpp FilePathParser.hpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

reindex_pdb: reindex_pdb.cpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

strip_sidechain: strip_sidechain.cpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

BackboneTorsion: BackboneTorsion.cpp ${HEADER}
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

clean:
	rm ${PROG}
