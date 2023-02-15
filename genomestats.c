#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFSIZE 1000000

int main(int argc, char **argv)
{
	FILE* iFPtr;
	
	if(argc != 2) {
		printf("Not enough args! Argc = %d.\n", argc);
		exit(1);
	}
	char* refSeqFile = argv[1];
	
	fprintf(stderr, "refSeqFile = %s.\n", refSeqFile);
	
	char* readBuffer = (char*)malloc(sizeof(char) * BUFSIZE);
	int len = 0;
	int lines = 0;
    unsigned long totalC = 0, totalG = 0;
	
	char* token;
	char seps[] = " \t\n\r";
	char chrName[100];
	
	if(!(iFPtr = fopen(refSeqFile, "r"))) {
		printf("File %s open error.\n", refSeqFile);
		exit(1);
	}
	
	while(fgets(readBuffer, BUFSIZE, iFPtr)) {
		if(strlen(readBuffer) >= BUFSIZE - 1) {
            fprintf(stderr, "Too many characters in one row! Try to split the long row into several short rows (fewer than %d characters per row).\n", BUFSIZE);
            exit(1);
        }
		
		if(readBuffer[0] == '>') {
			// Save name
			token = strtok(readBuffer + 1, seps);	
			strcpy(chrName, token);
			len = 0;	
		}
		else {
			// Substract \n
			//len += strlen(readBuffer) - 1;
            len=0;
            while(readBuffer[len]!='\0') {
                if(readBuffer[len]=='C' || readBuffer[len]=='c') totalC++;
                else if(readBuffer[len]=='G' || readBuffer[len]=='g') totalG++;
                len++;
            }
		}
		
		lines++;
	}
	
	fclose(iFPtr);
	
    fprintf(stdout, "C\tG\tCG\n%ld\t%ld\t%ld\n", totalC, totalG, totalC+totalG);
	return 0;
}
