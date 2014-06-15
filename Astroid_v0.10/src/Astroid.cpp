/*    
 *    Astroid.cpp		
 *    Astroid
 *
 *    Copyright (C) 2014 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define UNIX

#ifdef UNIX

#include <fstream>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>

#else

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <string>

#endif


using namespace std;

const int default_dataset_num = 100;


string itostr(long t);
void inputChrName(string namefile);
unordered_map<string, long> chrInfo;


int main(int argc, char *argv[])
{
	if (argc != 4)
	{
		cout << argv[0] << "<AstroidPath>\t<SAMFileName>\t<OutputFilePath>" << endl;
		return 1;
	}

	string comd, filepath, srcpath, filename;
	int i;
	srcpath = argv[1];


	/************************************************************************/
	// temporary file folder	
	filepath = argv[3];
	comd = "mkdir " + filepath + "tmp_foler";
	i = system(comd.c_str());


	/************************************************************************/
	// parse the sam file
	comd = "mkdir " + filepath + "tmp_foler/SAM";
	i = system(comd.c_str());

	filename = argv[2];
	comd = srcpath + "bin/sepSAM " + filename + " " + filepath + "tmp_foler/SAM/";
	i = system(comd.c_str());

	/************************************************************************/
	// sort the sam file
	filename = filepath + "tmp_foler/SAM/ChromosomeName.txt";
	inputChrName(filename);
	unordered_map<string, long> ::iterator iter = chrInfo.begin();
	while(iter != chrInfo.end())
	{
		comd = "sort -k 1,1n " + filepath + "tmp_foler/SAM/" + iter->first + ".txt > " + filepath + "tmp_foler/SAM/" + iter->first + "_sorted.txt";
		i = system(comd.c_str());
		comd = "rm -r " + filepath + "tmp_foler/SAM/" + iter->first + ".txt";
		i = system(comd.c_str());
		comd = "mv " + filepath + "tmp_foler/SAM/" + iter->first + "_sorted.txt " + filepath + "tmp_foler/SAM/" + iter->first + ".txt";
		i = system(comd.c_str());
		iter++;
	}

	cout << "Parse SAM file done..." << endl;

	/************************************************************************/
	// parse the fragments
	comd = "mkdir " + filepath + "tmp_foler/Frag";
	i = system(comd.c_str());
	iter = chrInfo.begin();
	while(iter != chrInfo.end())
	{
		comd = "mkdir " + filepath + "tmp_foler/Frag/" + iter->first;
		i = system(comd.c_str());
		comd = srcpath + "bin/fragment " + filepath + "tmp_foler/SAM/ " + iter->first + " 1 noDB " + filepath + "tmp_foler/Frag/" + iter->first + "/";
		i = system(comd.c_str());
		iter++;
	}

	cout << "Parse Fragments done..." << endl;

	/************************************************************************/
	// run Astroid
	comd = "mkdir " + filepath + "tmp_foler/Result";
	i = system(comd.c_str());
	iter = chrInfo.begin();
	while(iter != chrInfo.end())
	{
		comd = "mkdir " + filepath + "tmp_foler/Result/" + iter->first;
		i = system(comd.c_str());
		comd = "mkdir " + filepath + "tmp_foler/Result/" + iter->first + "/tmp";
		i = system(comd.c_str());
		comd = srcpath + "bin/Assembly " + iter->first + " " + filepath + "tmp_foler/Frag/" + iter->first + "/ " + filepath + "tmp_foler/SAM/ " + iter->first + " ";
		comd += filepath + "tmp_foler/Result/" + iter->first + "/tmp/ " + filepath + "tmp_foler/Result/" + iter->first + "/" + iter->first + ".gtf " + itostr(iter->second);
		i = system(comd.c_str());

		iter++;
	}
	

	/************************************************************************/
	// delete temporary results
	filename = filepath + "Astroid_result.txt";
	comd = "cat ";
	iter = chrInfo.begin();
	while(iter != chrInfo.end())
	{
		comd += filepath + "tmp_foler/Result/" + iter->first + "/" + iter->first + ".gtf ";
		iter++;
	}
	comd += "> " + filename;
	i = system(comd.c_str());
	
	comd = "rm -r " + filepath + "tmp_foler";
	i = system(comd.c_str());
	chrInfo.clear();

	return 0;
}

// Transform long to string
string itostr(long t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}

void inputChrName(string namefile)
{
	ifstream chrFile;
	string chrNm, info;
	long chrsize;
	chrFile.open(namefile.c_str());

	while(chrFile >> chrNm)
	{
		chrFile >> info;
		chrsize = atol(info.c_str());
		getline(chrFile, info);
		if (chrNm[0] == 'c' && chrNm[1] == 'h' && chrNm[2] == 'r')
		{
			chrInfo.insert(make_pair<string, long>(chrNm, chrsize));
		}

	}
	chrFile.close();

	return;
}