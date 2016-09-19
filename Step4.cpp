//generate Term Frequency . Remove one line in output after completion
#include<iostream>
#include<fstream>
#include<map>
using namespace std;
int main(){
	ifstream input,global;
	ofstream porterf;
	int count;
	int i=0;
	int num;
	string s;
	map<string,int> temp;
	map<string,int>globalmap; 
	string word;
	input.open("/home/livingroom/DC/output_files/Step1_final.txt");
	global.open("/home/livingroom/DC/output_files/global2.txt");	
	porterf.open("/home/livingroom/DC/output_files/TermFrequency.txt");
	//cout<<"---g"<<endl;
	while(global>>s){
		//cout<<i++<<endl;
			++globalmap[s];
	}
		//cout<<"---h"<<endl;
	 while(!input.eof()){
	 	input>>num;
	 	//cout<<num<<endl;
	 	if(num==0)
	 		porterf<<"0"<<endl;
	 	else if(num>0){
	 		while(num>0){
	 		input>>s;
	 		if(globalmap.find(s)!=globalmap.end())
	 			++temp[s];
	 			num--;
	 	}
	 	porterf<<temp.size()<<" ";
	 	//cout<<temp.size()<<" ";
	 	for (map<string, int>::iterator it = temp.begin(); it != temp.end(); ++it){	 					
	     	     	porterf<<it->first<<" "<<it->second<<" ";  		          	
	     	     	//cout<<it->first<<" ";  		          	
	    }
	    porterf<<endl;
	    temp.clear();
	    }
	 }
	 input.close();
	 global.close();
	 porterf.close();
	return 0;
}
