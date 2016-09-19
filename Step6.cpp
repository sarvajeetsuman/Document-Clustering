//finding inverse document frequency 
#include<iostream>
#include<cstdlib>
#include<map>
#include<cstring>
#include<fstream>
#include<cstdlib>
#include<vector>
#include <unistd.h>
#include<cmath>
using namespace std;

int main()
{
	int count=0;
	map<string,int>dictionary;
	map<string,int>wordCount;
	map<int,int>mat;
	std::map<std::string, int>::iterator it2;
	ifstream matfile;
	ofstream myfile;
	string line2,body,fullbody;
	int line;
	int p,q,r;
	float idf=0.0;
	myfile.open("/home/livingroom/DC/output_files/idf.txt");
	matfile.open("/home/livingroom/DC/output_files/mat.txt");
	
	while(matfile >> line){
		if(line<=0)
			continue;
		++mat[line];       	
	}
	for (map<int, int>::iterator it = mat.begin(); it != mat.end(); ++it){
	          idf=1+log(21578/mat.find(it->first)->second);
	        	myfile<<it->first<<" "<<idf<<endl;
	 }	
	myfile.close();
	matfile.close();
}
