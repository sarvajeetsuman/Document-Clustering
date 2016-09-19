// Snn & Random walk Based Model implemented in STXXL/c++ 
#include<iostream>
#include<limits>
#include<stxxl/vector>
#include<stxxl/stream>
#include<stxxl/bits/containers/matrix.h>
#include<fstream>
#include<string>


using namespace stxxl;
using namespace std;

typedef  double value_type;

int number_of_point=21578;
int threshold=7;
int k=0;


typedef stxxl::VECTOR_GENERATOR<unsigned int, 1, 1, 64*1024*1024, stxxl::RC, stxxl::lru>::result vector_type;
vector_type my_vector;


bool isInCorePoints(int element)
{
	vector_type::const_iterator iter = my_vector.begin();
		for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
			{
				if(*iter==(unsigned int)element)
					return true;
				//for(int z=0;z<k;z++)
				//{
					++iter;
					//j++;
				//}
				
			}
		return false;
}
void print()
{
	vector_type::const_iterator iter = my_vector.begin();
		for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
			{
				printf("%d",*iter);
				++iter;
				
			}
		
}


int main()
	{
		
		const int block_order =64;
		const int block_order2=64;
		ifstream doc;
		
		double idf;
		doc.open("/media/livingroom/stxxl/stxxl-1.4.1/local/myfile_final64.txt",ios::in); //read from input file
		
		char file[]="/media/livingroom/stxxl/stxxl-1.4.1/local/files.txt";//tell about document in which cluster 
		
	
		
		int size=21578; //number of document
		std::cout<<"size is "<<size<<"\n";
		
		
		int dimension=64; //dimension of document(no. of term)
		std::cout<<"dimension is "<<dimension<<"\n";
		
		k=64;//size of cluster
		std::cout<<"k is "<<k<<"\n";
		
		number_of_point=size;
		std::cout<<"block order is "<<block_order2<<"\n";
		int rows=floor(block_order2 * block_order2/dimension);
		std::cout<<"no of data item rows is "<<rows<<"\n";
		
		int rows_last=size%rows;
		int no_of_blocks=(size/rows);
		if(rows_last>0)
			no_of_blocks+=1;
		std::cout<<"no of blocks  is "<<no_of_blocks<<"\n";
		
		
		std::cout<<"no of rows in last block  is "<<rows_last<<"\n";
		
		typedef swappable_block_matrix<value_type, block_order> swappable_block_matrix_type;
		typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
		typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;
		typedef typename block_scheduler_type::internal_block_type internal_block_type;


	
		typedef swappable_block_matrix<value_type, block_order2> swappable_block_matrix_type2;
		typedef typename swappable_block_matrix_type2::block_scheduler_type block_scheduler_type2;
		typedef typename block_scheduler_type2::swappable_block_identifier_type swappable_block_identifier_type2;
		typedef typename block_scheduler_type2::internal_block_type internal_block_type2;
		
       
        typedef block_scheduler< matrix_swappable_block<value_type, block_order> > bst;
        typedef matrix<value_type, block_order> mt;
        
        typedef block_scheduler< matrix_swappable_block<value_type, block_order2> > bst2;
        typedef matrix<value_type, block_order2> mt2;
		
		bst * b_s=new bst(2*sizeof(value_type)*block_order*block_order);
		bst2 * b_s2=new bst2(2*sizeof(value_type)*block_order2*block_order2);
		bst &bs=*b_s;
		bst2 &bs2=*b_s2;
	
		mt2 * a = new mt2(bs2, no_of_blocks*block_order2,block_order2);
		mt * c = new mt(bs, no_of_blocks*block_order,block_order);
		
	
		STXXL_MSG("intializing Matrices\n");
		
		for(int row_block_no=0;row_block_no<no_of_blocks;row_block_no++)
			{	
		
				for(int col_block_no=0;col_block_no<1;col_block_no++)
					{
						
			
					//	printf("initializing %d and %d block of matrix A \n",row_block_no,col_block_no);
						internal_block_type2 &iba=((*((*a).data)).bs).acquire((*((*a).data))(row_block_no,col_block_no));
						int last_block=0;
						if(rows_last>0)
							last_block=no_of_blocks-1;
						else
							last_block=no_of_blocks;
						
						if(row_block_no!=last_block)
						{
							
			
			
							for(int i=0;i<rows;i++)
								{
									
									for(int j=0;j<dimension;j++)
										{
											doc>>idf; 
											iba[i*dimension + j]=idf;
													
										}
								
								
								}
			
						}
						else
						{
							
							for(int i=0;i<rows_last;i++)
								{
									for(int j=0;j<dimension;j++)
										{
											doc>>idf; 
											iba[i*dimension + j]=idf;
													
										}
								
								}	
						}
						((*((*a).data)).bs).release((*((*a).data))(row_block_no,col_block_no),true);
					}
		
			}
		
		
		printf("no_of_blocks=%d \n",no_of_blocks);
		for(int row_block_no=0;row_block_no<no_of_blocks;row_block_no++)
			{
			
				for(int col_block_no=0;col_block_no<1;col_block_no++)
					{
						
						internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(row_block_no,col_block_no));
						int last_block=0;
						if(rows_last>0)
							last_block=no_of_blocks-1;
						else
							last_block=no_of_blocks;
						if(row_block_no!=last_block)
						{
							
							for(int i=0;i<rows;i++)
								{
									for(int j=0;j<k;j++)
										{
											if(j==0)
												{
													
													ibD[i*dimension + j]=row_block_no *rows+i;
													
												}
											else
												{
													
													ibD[i*dimension + j]=999999;
													
												}
									
										
										}
										
								}
						}
						else
						{
							for(int i=0;i<rows_last;i++)
								{
									for(int j=0;j<=k;j++)
										{
											if(j==0)
												{
													ibD[i*dimension + j]=row_block_no *rows+i;
													
												}
											else
												{
													ibD[i*dimension + j]=999999;
													
												}
									
										
										}
										
								}
						}
						((*((*c).data)).bs).release((*((*c).data))(row_block_no,col_block_no),true);
					}
		
			}
			
				
		printf("Initialization is complete\n");
		printf("Number of points %d\n",number_of_point);
		printf("Number of Near neighbors %d\n",k);
        printf("Threshold is %d\n",threshold);
		
/************************************************ starting Proposed method*******************************************************/
		double sum=0;
		double p,pp;
		int flag=0,n=0,da;
		stats_data stats_before, stats_after;
	    matrix_operation_statistic_data matrix_stats_before, matrix_stats_after;
	    matrix_stats_before.set();
	    stats_before = *stats::get_instance();
	    printf("preparing the k-neighbourhood matrix\n");  
		int rowsTemp1,rowsTemp2;
		double **distanceMat = (double **) malloc(sizeof(double *)*rows);
		int i=0;//j=0;
		for(i=0; i<rows; i++)/* Allocate array, store pointer  */
			{
				distanceMat[i] = (double *) malloc(sizeof(double)*k); 
			}
			//out1<<"\n starting Proposed method\n";
		for(int acol=0;acol<1;acol++)
			{
				for(int arow=0;arow<no_of_blocks;arow++)
					{
						for(int kk=0;kk<rows;kk++)
							for(int ll=0;ll<k;ll++)
								distanceMat[kk][ll]=99999;
						
						internal_block_type2 &ibA=((*((*a).data)).bs).acquire((*((*a).data))(arow,acol));
						
						internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(arow,acol));
						if(arow==no_of_blocks-1)
							rowsTemp1=rows_last;
						else
							rowsTemp1=rows;
						for(int ar=0;ar<1;ar++)
							for(int brow=0;brow<no_of_blocks;brow++)
								{
									internal_block_type2 &ibB=((*((*a).data)).bs).acquire((*((*a).data))(brow,ar));
									if(brow==no_of_blocks-1)
										rowsTemp2=rows_last;
									else
										rowsTemp2=rows;
									for(int i=0;i<rowsTemp1;i++)
										{
											for(int t=0;t<rowsTemp2;t++)
												{
													double z=0;
													
													/*calculating distance*/
													for(int j=0;j<dimension;j++)
														{
															z = z+ ( ibA[i*dimension + j] - ibB[t*dimension + j] )*( ibA[i*dimension + j] - ibB[t*dimension + j] );
														}
													sum = sqrt(z);
													
													if(distanceMat[i][k-1]> sum)
														{
															
															distanceMat[i][k-1] = sum;
													
															ibD[i*(k+1) + k] = t+brow*rows;	  //same problem
															da=k-1;
															while ( da >0 &&  distanceMat[i][da] < distanceMat[i][da-1]) 
																{
																	p=distanceMat[i][da];
																	pp= ibD[i*(k+1)+ da+1];
																	distanceMat[i][da]=distanceMat[i][da-1];
																	ibD[i*(k+1)+ da+1]   = ibD[i*(k+1) + da];
																	distanceMat[i][da-1]=p;
																	ibD[i*(k+1) + da] = pp;
	 																
																	da--;
																}
														}
												}
										}
									((*((*a).data)).bs).release((*((*a).data))(brow,ar),false);
								}
	
							
							((*((*a).data)).bs).release((*((*a).data))(arow,acol),false);
							((*((*c).data)).bs).release((*((*c).data))(arow,acol),true);
					}
				
	
	
			}
			
	
	
		//for(int row_block_no=0;row_block_no<no_of_blocks;row_block_no++)
			//{
			
				//for(int col_block_no=0;col_block_no<1;col_block_no++)
					//{
				
						//internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(row_block_no,col_block_no));
						//if(row_block_no!=(no_of_blocks-1))
						//{
							//for(int i=0;i<rows;i++)
								//{
									//for(int j=0;j<=k;j++)
										//{
										
													//printf("%lf",ibD[i*(k+1) + j]);
												
									
										
										//}
										//printf("\n");
								//}
						//}
						//else
						//{
							//for(int i=0;i<rows_last;i++)
								//{
									//for(int j=0;j<=k;j++)
										//{
										
													//printf("%lf",ibD[i*(k+1) + j]);
												
									
										
										//}
										//printf("\n");
								//}
						//}
						//((*((*c).data)).bs).release((*((*c).data))(row_block_no,col_block_no),true);
					//}
		
			//}
	printf("starting  Proposed method\n");
	
	stxxl::stats_data begin_1(*stats::get_instance());	
		
		/********************************************
		  calculating number of shared near neighbors.
		********************************************/
		
		int minPts=10;
		for(int acol=0;acol<1;acol++)
			{ 	
				int *sNN=(int *)malloc(sizeof(int)*rows);
				int *changed=(int *)malloc(sizeof(int)*rows);
				for(int arow=0;arow<no_of_blocks;arow++)
					{	
					
						for(int xx=0;xx<rows;xx++)
						{
							sNN[xx]=0;
							changed[xx]=0;
						}
						internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(arow,acol));
						if(arow==no_of_blocks-1)
							rowsTemp1=rows_last;
						else
							rowsTemp1=rows;
						for(int col=0;col<1;col++)
							{
								for(int row=0;row<no_of_blocks;row++)
									{
										internal_block_type &ibd=((*((*c).data)).bs).acquire((*((*c).data))(row,col));
										if(row==no_of_blocks-1)
											rowsTemp2=rows_last;
										else
											rowsTemp2=rows;
										
										for(int i=0;i<rowsTemp1;i++)
											{
												for(int t=0;t<rowsTemp2;t++)
													{
														n=0;
														
														for(int j=1;j<=k &&flag!=1;j++)
															{
																if(ibD[i*(k+1) + 1]==ibd[t*(k+1)+ j])
																	flag=1;
																
															}
														for(int j=1;j<=k &&flag!=2;j++)
															{
																if(ibD[i*(k+1) + j]==ibd[t*(k+1) + 1])
																	flag=flag+1;
															}
														if(flag==2)
															{
																for(int rr=1;rr<=k;rr++)
																	for(int yy=1;yy<=k;yy++)
																		{
																			if(ibD[i*(k+1) + rr]==ibd[t*(k+1) + yy])
																				n=n+1;
																		}
															}
	
														
														if(n>threshold)
															{// printf("here %d\t",rows*row+t);
															//printf("size  %d\t",my_vector.size());
																sNN[i]=sNN[i]+1;
																if(isInCorePoints(row*(rows)+t))
																{
																	
																	if(changed[i]==0)
																	{
																		
																		ibD[i*(k+1)+0]=row*rows+t;
																		changed[i]=1;
																		}
																}
	
															}
														
														n=0;
														flag=0;
													}
											}
										((*((*c).data)).bs).release((*((*c).data))(row,col),false);
									}
							}
			
						if(arow==no_of_blocks-1)
							rowsTemp1=rows_last;
						else
							rowsTemp1=rows;
						for(int itr=0;itr<rowsTemp1;itr++)
							
								{
									if(sNN[itr]>=minPts)
										{
											//for(int x=0;x<=k;x++)
											//{
												my_vector.push_back(arow*rows+itr);									
												
											//}
										}
									else 
										{
											if(changed[itr]==1)											
												ibD[itr*(k+1)+0]=arow*rows+itr;										
										}
								
								}
						
							((*((*c).data)).bs).release((*((*c).data))(arow,acol),true);	
					}
				
			}
			
			
		printf("no of core Points : %u\n",my_vector.size());
	
	//for(int row_block_no=0;row_block_no<no_of_blocks;row_block_no++)
			//{
			
				//for(int col_block_no=0;col_block_no<1;col_block_no++)
					//{
				
						//internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(row_block_no,col_block_no));
						//if(row_block_no!=(no_of_blocks-1))
						//{
							//for(int i=0;i<rows;i++)
								//{
									//for(int j=0;j<=k;j++)
										//{
										
													//printf("%lf\t",ibD[i*(k+1) + j]);
												
									
										
										//}
										//printf("\n");
								//}
						//}
						//else
						//{
							//for(int i=0;i<rows_last;i++)
								//{
									//for(int j=0;j<=k;j++)
										//{
										
													//printf("%lf\t",ibD[i*(k+1) + j]);
												
									
										
										//}
										//printf("\n");
								//}
						//}
						//((*((*c).data)).bs).release((*((*c).data))(row_block_no,col_block_no),true);
					//}
		
			//}
	  
	   
		printf("Now assigning points to nearest core point\n");
	/*****************************************************************************
	 * Code for assigning points to any core points decided in previous step.
	 * **************************************************************************/
	 int index=0,temp=0;
	 vector_type::const_iterator tmp;
	 for(int acol=0;acol<1;acol++)
			{ 	
				int *sNN=(int *)malloc(sizeof(int)*rows);
				
				for(int arow=0;arow<no_of_blocks;arow++)
					{	
						//std::fill(sNN,sNN+sizeof(sNN),0);
						for(int xx=0;xx<rows;xx++)
							sNN[xx]=0;
						int flag=-1;
						if(arow==no_of_blocks-1)
							rowsTemp1=rows_last;
						else
							rowsTemp1=rows;
						internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(arow,acol));
						
					    vector_type::const_iterator iter = my_vector.begin();
						for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
							{
								
								
								int corePoint=*iter;
								iter++;
							
								index=corePoint/rows ;
								temp=corePoint%(rows);
							
								
								internal_block_type &ibd=((*((*c).data)).bs).acquire((*((*c).data))(index,acol));
								
								//tmp=iter;						
								for(int i=0;i<rowsTemp1;i++)
									{
									//	print();
										//if(!isInCorePoints(arow*rows+i))
											
										{	
										
										
										n=0;flag=0;
										for(int j=1;j<=k && flag!=1;j++)
											{
												
												if(ibD[i*(k+1) + 1]==ibd[temp*(k+1) + j])
													flag=1;
													
											}
										for(int j=1;j<=k &&flag!=2;j++)
											{
												if(ibD[i*(k+1) + j]==ibd[temp*(k+1) + 1])
													flag=flag+1;
											}
										if(flag==2)
											{
												for(int rr=1;rr<=k;rr++)
													for(int yy=1;yy<=k;yy++)
														{
															if(ibD[i*(k+1) + rr]==ibd[temp*(k+1) + yy])
																n=n+1;
														}
											}
										if (n>=threshold)
											{
												if (n>sNN[i])
													{ 
														ibD[i*(k+1)+0]=ibd[temp*(k+1)+0];
														sNN[i]=n;
													}
											}
										
										
										
									}}
									((*((*c).data)).bs).release((*((*c).data))(index,acol),false);
							}
							
							((*((*c).data)).bs).release((*((*c).data))(arow,acol),true);
					}	
			}
	 /***
	  * 
	  * ********/
	  printf("This is second");
	  		stxxl::stats_data begin_2(*stats::get_instance());
	          	
	   STXXL_MSG(begin_2-begin_1);
	
		stats_after = *stats::get_instance();
	    matrix_stats_after.set();
	    STXXL_MSG(stats_after - stats_before);
	    
	    //for(int row_block_no=0;row_block_no<no_of_blocks;row_block_no++)
			//{
			
				//for(int col_block_no=0;col_block_no<1;col_block_no++)
					//{
				
						//internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(row_block_no,col_block_no));
						//if(row_block_no!=(no_of_blocks-1))
						//{
							//for(int i=0;i<rows;i++)
								//{
									//for(int j=0;j<=k;j++)
										//{
										
													//printf("%lf",ibD[i*(k+1) + j]);
												
									
										
										//}
										//printf("\n");
								//}
						//}
						//else
						//{
							//for(int i=0;i<rows_last;i++)
								//{
									//for(int j=0;j<=k;j++)
										//{
										
													//printf("%lf",ibD[i*(k+1) + j]);
												
									
										
										//}
										//printf("\n");
								//}
						//}
						//((*((*c).data)).bs).release((*((*c).data))(row_block_no,col_block_no),true);
					//}
		
			//}
	  
	  int count=0;
	  for(int acol=0;acol<1;acol++)
	{
		for(int arow=0;arow<no_of_blocks;arow++)
		{
			if(arow==no_of_blocks-1)
				rowsTemp1=rows_last;
			else
				rowsTemp1=rows;
			internal_block_type &ibc=((*((*c).data)).bs).acquire((*((*c).data))(arow,acol));
	
				for(int i=0;i<rowsTemp1;i++)
				{ 
					if(ibc[i*(k+1)+0]!=ibc[i*(k+1)+1])
					{
					//	printf("%f\n",ibc[i*element_row+1]);
						count++;
					}
				}
		((*((*c).data)).bs).release((*((*c).data))(arow,acol),false);
		}
	}
	
	int num_clusters=0;
	vector_type::const_iterator iter = my_vector.begin();
						for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
							{
								
								int corePoint=*iter;
								iter++;
								
								index=corePoint/rows;
								temp=corePoint%(rows);
															
								internal_block_type &ibd=((*((*c).data)).bs).acquire((*((*c).data))(index,0));
								
								if(ibd[temp*(k+1) + 1]==ibd[temp*(k+1)+ 0]){
									
									num_clusters++;
								}
								//cout<<ibd[temp*(k+1) + 1] <<"="<<ibd[temp*(k+1)+ 0]<<"NumClusters ="<<num_clusters<<endl;
								
									((*((*c).data)).bs).release((*((*c).data))(index,0),false);
							
	
								
								
							}
	
	printf("number of clusters: %d\n",num_clusters);
	printf("no of points assigned :%d\n",count);
	printf("no of core points %d\n",my_vector.size());
	
	std::ofstream labelFile;
	
	labelFile.open(file);
	
	for(int row_block_no=0;row_block_no<no_of_blocks;row_block_no++)
			{
			
				for(int col_block_no=0;col_block_no<1;col_block_no++)
					{
				
						internal_block_type &ibD=((*((*c).data)).bs).acquire((*((*c).data))(row_block_no,col_block_no));
						if(row_block_no!=(no_of_blocks-1))
						{
							for(int i=0;i<rows;i++)
								{
									labelFile<<ibD[i*(k+1)+0]<<"\n";
										
									//printf("%lf",ibD[i*(k+1) + 0]);
									//printf("\n");
								}
						}
						else
						{
							for(int i=0;i<rows_last;i++)
								{
									labelFile<<ibD[i*(k+1)+0]<<"\n";
									//printf("%lf",ibD[i*(k+1) + 0]);
									//printf("\n");
								}
						}
						((*((*c).data)).bs).release((*((*c).data))(row_block_no,col_block_no),true);
					}
		
			}
			labelFile.close();
	  
	/**********************
	 * 
	 * ***********/
		

		
	/********************************************** END OF THE PROGRAM **********************************************************/
		
		
	
		STXXL_MSG("END OF THE PROGRAM\n");
		
		my_vector.clear();
		doc.close();
		
		return 1;
	}

