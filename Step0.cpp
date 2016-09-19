//this file has been created to collect the text & body in 1 file
#include<iostream>
#include<ctype.h>
#include<vector>
#include<cstring>
#include<map>
#include<cstdlib>
#include<fstream>
using namespace std;
int main()
{
	 int i;
     ifstream input;
     ofstream myfile;
     bool flag;
     int rcount=0;
     string line,date,title,body,fullbody;
	 int p,q,r;
  	 myfile.open ("/home/livingroom/DC/output_files/Step0.txt",ios::app);
  	 
     for(i=0;i<22;i++){
	 	char base[100]="/home/livingroom/DC/reuters21578/reut2-";
   		char extension[5]=".sgm";
	   	char filenumber[5];
		sprintf(filenumber,"%03d",(i));
	    strcat(base,strcat(filenumber,extension));
	    cout<<base<<endl;
	     
	     input.open(base);
	     
	     while(input.good())
	     {
	     	getline(input,line);
	     	flag=false;
			 q=line.find("<REUTERS");
			 if(q>=0){
			 rcount++;
			 	myfile<<"\n"<<rcount<<" ";
			 	fullbody="";
			 	line=line.substr(q);
			 	/*while(line.find("<TITLE>")==-1 && line.find("</REUTERS>")==-1)
			 	{
			 		getline(input,line);
			 	
			 	}
			 	if(line.find("</REUTERS>")==-1)
			 	{
			 		myfile<<line<<endl;
			 	}
			 	else{
			 		myfile<" NULL ";
			 	}*/
			 	while(line.find("<BODY>")==-1 && line.find("</REUTERS>")==-1)
			 	{
			 		getline(input,line);
			 	
			 	}
			 	if(line.find("<BODY>")>=0 && line.find("</REUTERS>")==-1)
			 	{
			 		flag=true;
			 	p=line.find("<BODY>");
			 	line=line.substr (p);
			 	//body=line.replace(line.find("<BODY>"),6,"");
			 	r=line.find("</BODY>");
			 	while(r==-1){
			 		myfile<<line<<endl;
			 		//cout<<line;
			 		fullbody.append(line);
			 		getline(input,line);
			 		r=line.find("&#3;</BODY></TEXT>");	
			 	}
			 	myfile<<"</BODY>"<<endl;
			 	}
			 	else{
			 		myfile<<"<BODY> NULLNULL </BODY>"<<endl;
			 	}
	     }
	    
	}
	input.close();
}
myfile.close();
}

