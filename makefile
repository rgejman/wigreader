all: FCWig wigToBed geneProfile profileFromEnds genomicColocalization subtractWig promoterColocalization normalizeWig intergenicProfile

FCWig: wigreader.o FCWig.o functions.o
	g++ -o bin/FCWig wigreader.o FCWig.o functions.o -lpthread 
	
FCWig.o: FCWig.cpp wigreader.h functions.h
	g++ -c FCWig.cpp

subtractWig: wigreader.o subtractWig.o functions.o
	g++ -o bin/subtractWig wigreader.o subtractWig.o functions.o -lpthread 
	
	
subtractWig.o: subtractWig.cpp wigreader.h functions.h
	g++ -c subtractWig.cpp


wigToBed: wigreader.o wigToBed.o functions.o
	g++ -o bin/wigToBed wigreader.o wigToBed.o functions.o
	
	
wigToBed.o: wigToBed.cpp wigreader.h functions.h
	g++ -c wigToBed.cpp

intergenicProfile: wigreader.o intergenicProfile.o functions.o
	g++ -o bin/intergenicProfile wigreader.o intergenicProfile.o functions.o
	
	
intergenicProfile.o: intergenicProfile.cpp wigreader.h functions.h
	g++ -c intergenicProfile.cpp

normalizeWig: wigreader.o normalizeWig.o functions.o
	g++ -o bin/normalizeWig wigreader.o normalizeWig.o functions.o
	
	
normalizeWig.o: normalizeWig.cpp wigreader.h functions.h
	g++ -c normalizeWig.cpp

promoterColocalization: wigreader.o promoterColocalization.o functions.o
	g++ -o bin/promoterColocalization wigreader.o promoterColocalization.o functions.o -lpthread 
	
promoterColocalization.o: promoterColocalization.cpp wigreader.h functions.h
	g++ -c promoterColocalization.cpp

genomicColocalization: wigreader.o genomicColocalization.o functions.o
	g++ -o bin/genomicColocalization wigreader.o genomicColocalization.o functions.o -lpthread 
	
genomicColocalization.o: genomicColocalization.cpp wigreader.h functions.h
	g++ -c genomicColocalization.cpp

geneProfile: wigreader.o geneProfile.o functions.o
	g++ -o bin/geneProfile wigreader.o geneProfile.o functions.o
	
geneProfile.o: geneProfile.cpp wigreader.h functions.h
	g++ -c geneProfile.cpp
	
profileFromEnds: wigreader.o profileFromEnds.o functions.o
	g++ -o bin/profileFromEnds wigreader.o profileFromEnds.o functions.o -lpthread 
	
profileFromEnds.o: profileFromEnds.cpp wigreader.h functions.h
	g++ -c profileFromEnds.cpp

wigreader.o: wigreader.cpp wigreader.h functions.h
	g++ -c wigreader.cpp
	
functions: functions.o
	g++ -o bin/functions functions.o
	
functions.o: functions.h
	g++ -c functions.cpp
	
clean:
	rm bin/FCWig bin/profileFromEnds bin/geneProfile bin/genomicColocalization bin/promoterColocalization bin/normalizeWig bin/intergenicProfile bin/subtractWig *.o