/*
This function is an aid to compiling PARI libraries from separate files.
Assume you have three files f1.c, f2.c, and f3.c that you want to compile into one library, named "allfcns".
1.	Compile this function with "gcc pari_compile.c -o pari_compile"
2.	Create a text file (say makefile.txt) which contains:
		3
		allfcns
		f1 f2 f3
	The first line is the number of files, the second line is the library name, and the third line lists all files.
3.	Call "./pari_compile makefile" to compile the .o files (will be named f1.o, f2.o, and f3.o), and then make the shared library. The shared library is named liballfcns.so.

If you already have some of the .o files and only want to recompile some of them, add them to step 3. For example, "./pari_compile makefile f2 f3" will not touch f1.o, compile f2.o and f3.o, and create the shared library liballfcns.so from f1 f2 and f3. When installing functions in gp, the fourth entry (the library) will be "./liballfncs.so"
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//The only method: make the .o files, compile the library
void main(int nargs, char *args[]){
  int maxfilelen=50;
  if(nargs==1){
	printf("Pass as arguments: the file you store the file names to be compiled, and a list of all object files you want recompiled. Pass 0 for the second entry if you want ALL to be recompiled.\n");
	return;
  }
  char data[strlen(args[1])+4];
  sprintf(data, "%s.txt", args[1]);
  FILE *f=fopen(data, "r");//Open the file to get data from.
  int nfiles, s;
  char libraryname[maxfilelen];//The name of the library
  fscanf(f, "%d %s", &nfiles, libraryname);//The number of files, and the name of the library
  char libcreate[143+(maxfilelen+3)*nfiles];//Allocating space, accounting for the space and the .o
  sprintf(libcreate, "gcc -o lib%s.so -shared -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC -Wl,-shared -Wl,-rpath=\"./\" -lc -lm -L/usr/local/lib -lpari", libraryname);//The final call to create the library
  char objectcommand[98+2*maxfilelen];
  char eachfile[50];//Filenames are 50 characters or less
  if(nargs>2){//We specify the .o files created
    for(int i=2;i<nargs;i++){//Creating the .o files
	  sprintf(objectcommand, "gcc -c -o %s.o -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC -I\\\"/usr/local/include\\\" %s.c", args[i], args[i]);
	  s=system(objectcommand);
	  if(s==-1){printf("ERROR creating %s object file\n", args[i]);fclose(f);return;}
	}
	for(int i=0;i<nfiles;i++){//For each file, we append it to the libcreate
      fscanf(f, "%s", eachfile);//Read the file
	  strcat(libcreate, " ");
	  strcat(libcreate, eachfile);
	  strcat(libcreate, ".o");
    }  
  }
  else{//We create ALL object files
    for(int i=0;i<nfiles;i++){//For each file, we create the object file
      fscanf(f, "%s", eachfile);//Read the file
	  sprintf(objectcommand, "gcc -c -o %s.o -O3 -Wall -fno-strict-aliasing -fomit-frame-pointer -fPIC -I\\\"/usr/local/include\\\" %s.c", eachfile, eachfile);
	  strcat(libcreate, " ");
	  strcat(libcreate, eachfile);
	  strcat(libcreate, ".o");
	  s=system(objectcommand);
	  if(s==-1){printf("ERROR creating %s object file\n", eachfile);fclose(f);return;}
    }
  }
  s=system(libcreate);
  if(s==-1){printf("ERROR creating the library\n");fclose(f);return;}
  fclose(f);
  return;
}
