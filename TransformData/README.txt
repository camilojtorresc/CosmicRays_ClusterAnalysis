
Corsika Reader
==============

	This file is created to read the Corsika binary data and transform it to
	.txt files to do the cluster analysis without using ROOT 5.


Prepare your Corsika environment
================================

	To transform your binary data file you need to use the coast library. To
	active coast library you need to do the following steps:

	1. cd 'PATH-TO-CORSIKA-DIR'/corsika-77500
	2. ./coconut


	when you run ./coconut you install a Corsika executable to do the simul-
	ations with your preferred model options. To activate the coast library
	when coconut says:
	
	options:   'YOUR-MODEL-OPTIONS'

	Which additional CORSIKA program options do you need ?
	

	choose the “d1” option:
	"d1 - Inclined observation plane"

	and finish your installation. Now, ou are able to use the coast library.

Prepare your Makefile executable
================================

	To transform your binary file you need to have the following files in
	the same directory:

	CorsikaReader.cc
	CorsikaReader.o
	Makefile
	Binary_File

	Note: You need to change a code line inside the MakeFile as follows:

	ifndef COAST_DIR
	  COAST_DIR='PATH-TO-CORSIKA-DIR'/corsika-77500
	endif


	Now you are ready to transform your binary file.

Transform your Binary_File
==========================

	To transform your binary file you only need to run the next command line

	1. make reader


	That creates an executable file CorsikaReader  to run and transform your
	Binary_File.
	
	Then execute:

	2. ./CorsikaReader Binary_File


	This output two .txt files with the particle and showers information

	Binary_File_particle.txt
	Binary_File_showers.txt

