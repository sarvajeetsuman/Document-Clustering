// Frequency and matrix Generation for IDF 
#include<iostream>
#include<map>
#include<fstream>
using namespace std;
int main()
{
	ifstream step1f,global;
	ofstream mat,freq;
	int count=0;
	int num;
	string s;
	map<string,int> globalmap;
	step1f.open("/home/livingroom/DC/output_files/Step1_final.txt");
	global.open("/home/livingroom/DC/output_files/global2.txt");
	mat.open("/home/livingroom/DC/output_files/mat.txt");
	freq.open("/home/livingroom/DC/output_files/freq.txt");
	while(global>>s){		
			count++;
			globalmap[s]=count;
	}
	
	while(!step1f.eof()){
		step1f>>num;
		freq<<num;
		freq<<endl;
		if(num==0)
	 		mat<<"-1";
	 	else if(num>0){
		while(num>0){
			step1f>>s;
			if(globalmap.find(s)!=globalmap.end()){
				mat<<globalmap.find(s)->second<<" ";
				//cout<<globalmap.find(s)->second;
			}
			--num;
			//cout<<"num= "<<num<<endl;
		}
		} 
		mat<<endl;
	}
	step1f.close();
	mat.close();
	global.close();	
}
