//finding TFIDF 
#include<iostream>
#include<map>
#include<cstdlib>
#include<cstring>
#include<fstream>
#include<cstring>
#include<vector>
#include <unistd.h>
using namespace std;
int main(){
	ifstream global,idf,tf;
	ofstream tfidfcomp,numofterms,idf_file;
	int num , fre;
	map<string,int>dictMap;
	map<int,float>idfMap;
	int i=0,count,index;
	int id;
	float idfNum,IDF,TF;
	
	string w,word,line;  
	idf.open("/home/livingroom/DC/output_files/idf.txt");
	global.open("/home/livingroom/DC/output_files/global2.txt");
	tf.open("/home/livingroom/DC/output_files/TermFrequency.txt");
	numofterms.open("/home/livingroom/DC/output_files/numofterms.txt");
	
	count=0;
	while(global>>line){
		count++;  
		dictMap[line]=count;
	}
	 
	while(idf >> id){
		if(idf >> idfNum)
			idfMap[id]= idfNum;
	}
	/*
	for (map<string, int>::iterator it = dictMap.begin(); it != dictMap.end(); ++it){

	          	index=dictMap.find(it->first)->second; 
	          	cout<<"index"<<"="<<index<<endl;
	          	if(idfMap.find(index)!=idfMap.end())
	  				IDF=idfMap.find(index)->second;
	  			cout<<"idf"<<"="<<IDF<<endl;
	   }*/
	for(i=1;i<=21578;i++){
		char base[100]="/home/livingroom/DC/tfidffiles/tfidf-";
   	 	char extension[5]=".dat";
	 	char filenumber[5];
		sprintf(filenumber,"%05d",(i));
		strcat(base,strcat(filenumber,extension));
		cout<<base<<endl;
		idf_file.open(base);
		tf >> num;
		
		numofterms<<num<<" ";
		for(int j=0;j<num;j++){
			tf>> w;
			tf>> fre;
			TF=(float)fre/num;
			index=dictMap.find(w)->second;
			if(idfMap.find(index)!=idfMap.end())
				IDF=idfMap.find(index)->second;
			//cout<<"TF "<<TF<<endl;
			//cout<<"IDF "<<IDF<<endl;
			if(TF*IDF>0){
				//cout<<"TF IDF  "<<TF<<" "<<IDF<<endl;
				idf_file<<dictMap.find(w)->second<<" "<<TF*IDF<<" "; 
			}
		}
			idf_file.close();		
		
	}
	numofterms.close();

	global.close();
	tfidfcomp.close();
}
