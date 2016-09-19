//Finding Global Frequency of words after porter
#include<iostream>
#include<fstream>
#include<map>
using namespace std;
int main(){
	ifstream input;
	ofstream global;
	int count;
	int max=0;
	map<string,int> globalmap;
	string word;
	input.open("/home/livingroom/DC/output_files/Step1_final1.txt");	
	global.open("/home/livingroom/DC/output_files/global2.txt");
	
	 while(input >>word){
	 	++globalmap[word];
	 }
	 for (map<string, int>::iterator it = globalmap.begin(); it != globalmap.end(); ++it){	 			
	 			if(it->second>max)
	 				max=it->second;	 			
	 			if(it->second>2)
	     	     	global<<it->first<<" "<<endl;//<<it->second<<endl;  		          	
	        }
	        cout<<"Max = "<<max<<endl;
	return 0;
}
