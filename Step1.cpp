//this file has been made for removing punctuation marks & stopwords

#include<iostream>
#include<ctype.h>
#include<vector>
#include<cstring>
#include<map>
#include<cstdlib>
#include<fstream>
#include<string>
  
#include <utility>
#include <stdio.h>
#include <locale>
using namespace std;


std::string toLower(const std::string& s)
{
    std::string result;

    std::locale loc;
    for (unsigned int i = 0; i < s.length(); ++i)
    {
        result += std::tolower(s.at(i), loc);
    }
   
    return result;
}
int isalphabet(string s){
for(int i=0;i<s[i]!='\0';++i){
	if(!isalpha(s[i]))
		return 0;
		
		}
	return 1;
}

struct letter_only: std::ctype<char> 
{
    letter_only(): std::ctype<char>(get_table()) {}

    static std::ctype_base::mask const* get_table()
    {
        static std::vector<std::ctype_base::mask> 
            rc(std::ctype<char>::table_size,std::ctype_base::space);

        std::fill(&rc['A'], &rc['z'+1], std::ctype_base::alpha);
        return &rc[0];
    }
};

int main()
{
	 int i=0;
	 int count=0;
     std::map<std::string, int> wordCount;
     map<string,int> hs;
	 std::locale loc;
     ifstream stopword;
     std::ifstream input,myfile2;
     std::ofstream myfile1,tfmatrix;
     string line,body,fullbody;
	 int p,q,r;
	 std::string word;
	 stopword.open("/home/livingroom/DCFinal/stopword.txt");
	 while(stopword >> word)
	     {
	   		word=toLower(word);
	     	 ++hs[word];
	         //cout<<word<<endl;
	     }
	 input.imbue(std::locale(std::locale(), new letter_only())); //enable reading only letters!
	 input.open("/home/livingroom/DC/output_files/Step0.txt");
	 tfmatrix.open("/home/livingroom/DC/output_files/Step1.txt",ios::app);
	
	 while(input.good()){
	     getline(input,line);
		 q=line.find("<BODY>");
		 if(q>=0){
			 myfile1.open ("/home/livingroom/DC/temp/temp.txt");
			 fullbody="";
			 line=line.substr(q);
			 body=line.replace(line.find("<BODY>"),6,"");
			 r=line.find("</BODY>");
			 while(r==-1){
			 	
			 	myfile1<<line<<endl;
			 	fullbody.append(line);
			 	getline(input,line);
			 	r=line.find("</BODY>");	
			 	
			 }
			 myfile1.close();
			 myfile2.open ("/home/livingroom/DC/temp/temp.txt");
	     	 while(myfile2 >> word){
	     	
				word=toLower(word);
	     	 	if(hs.find(word)==hs.end()&& isalphabet(word) ==1 && word.length()>=5 ){
	        		++wordCount[word];
	        	}
	     	}
	     	tfmatrix<<wordCount.size()<<" ";
	     	count=0;
	     	for (std::map<std::string, int>::iterator it = wordCount.begin(); it != wordCount.end(); ++it){
	           	tfmatrix<< it->first <<" ";
	          	count++;
	        }
	        tfmatrix<<endl;
	        wordCount.clear();
	        myfile2.close();
	     	}
	     }
	     tfmatrix.close();
}
