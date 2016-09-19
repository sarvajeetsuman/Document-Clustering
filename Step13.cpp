//Kmeans Algorithm used for Comparison
#include<limits>
#include<stxxl/vector>
#include<stxxl/stream>
#include<stxxl/bits/containers/matrix.h>
#include<fstream>
#include<string>

#include <set>
#include <vector>
#include <iostream>

#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>


#include <stdlib.h>
#include<math.h>


using namespace stxxl;
using namespace std;
typedef double value_type ;


stxxl::unsigned_type block =64;
stxxl::unsigned_type threshold=7,minPts=10;


typedef stxxl::VECTOR_GENERATOR<unsigned int, 1, 1, 64*1024*1024, stxxl::RC, stxxl::lru>::result vector_type;
vector_type my_vector;



bool isInCorePoints(int element)
	{
		vector_type::const_iterator iter = my_vector.begin();
		for (stxxl::uint64 j = 0; j < my_vector.size(); j++)
			{
				if(*(iter)==element)
					return true;
				//for(int z=0;z<k;z++)
				//{
					++iter;
					//j++;
				//}
				
			}
		return false;
	}
int similarity(double v1[],double v2[],int size )
	{ 
		
		int flag=0,x=0,y=0,n=0;
		for( x=0;x<size && flag!=1;x++)
			{
				if(v1[0]==v2[x])
					flag=1;
			}
		for( x=0;x<size && flag!=2;x++)
			{
				if(v2[0]==v1[x])
					flag=flag+1;
			}
		if(flag==2)
			{
				for(x=0;x<size;x++)
					for(y=0;y<size;y++)
					{
						if(v1[x]==v2[y])
							n++;
						
					}
			}
			return n;
	}
		
double distance2(double a[],double b[],int width)
	{
		double sum=0;
		for(int i=0;i<width;i++)
		{
			sum+=(a[i]-b[i])*(a[i]-b[i]);
		}
		return sqrt(sum);
	}



		char file[]="/media/livingroom/stxxl/stxxl-1.4.1/local/files.txt";
		ifstream doc;
		double idf;
		
		
		int size=21578;
		
		int dimension=64;
		
		int k=64;
		
		int init=size*0.1;
		
		
		const int large_block_order=64;
		stxxl::unsigned_type no_of_vector=size;
		int rows=floor(large_block_order*large_block_order/dimension);
	
		const int small_block_order=64;
		int no_of_blocks=ceil(size/rows);

		
		stats_data stats_begin, stats_end, stats_cumulative;
	
		int internal_memory1 = 2 * sizeof(value_type)*large_block_order*large_block_order;
		int internal_memory2 = 2* sizeof(value_type)*small_block_order*small_block_order;
	
		typedef block_scheduler< matrix_swappable_block<value_type, small_block_order> > bst;
        typedef matrix<value_type, small_block_order> mt;
        
        typedef block_scheduler< matrix_swappable_block<value_type, large_block_order> > bst2;
        typedef matrix<value_type, large_block_order> mt2;
	 	
		bst * b_s=new bst(internal_memory2);
		bst2 * b_s2=new bst2(internal_memory1);
		bst &bs=*b_s;
		bst2 &bs2=*b_s2;
	
	// Create 2 matrices with given dimensions
		mt2 *a=new mt2(bs2, no_of_vector, dimension);
		mt *b=new mt(bs,no_of_vector,k+1);
		
		
		typedef mt2::row_major_iterator row_iterator;
		typedef mt::row_major_iterator row_iterator2;
		
		
	    
		
namespace Clustering{

  typedef double Coord;            // a coordinate
  typedef double Distance;         // distance
  typedef unsigned int Dimensions; // how many dimensions
  typedef unsigned int PointId;    // the id of this point
  typedef unsigned int ClusterId;  // the id of this cluster

  typedef std::vector<Coord> Point;    // a point (a centroid)
  typedef std::vector<Point> Points;   // collection of points

  typedef std::set<PointId> SetPoints; // set of points
  
  // ClusterId -> (PointId, PointId, PointId, .... )
  typedef std::vector<SetPoints> ClustersToPoints;
  // PointId -> ClusterId
  typedef std::vector<ClusterId> PointsToClusters; 
  // coll of centroids
  typedef std::vector<Point> Centroids;

  //
  // Dump a point
  //
  std::ostream& operator << (std::ostream& os, Point& p);

  //
  // distance between two points
  //
  Distance distance(const Point & x, const Point & y);

  //
  // Dump collection of Points
  //
  std::ostream& operator << (std::ostream& os, Points& cps);

  //
  // Dump a Set of points
  //
  std::ostream& operator << (std::ostream& os, SetPoints & sp);
    
  //
  // Dump centroids
  //
  std::ostream& operator << (std::ostream& os, Centroids & cp);
  

  //
  // Dump ClustersToPoints
  //
  std::ostream& operator << (std::ostream& os, ClustersToPoints & cp);

  //
  // Dump ClustersToPoints
  //
  std::ostream& operator << (std::ostream& os, PointsToClusters & pc);


  //
  // This class stores all the points available in the model
  //
  class PointsSpace{

    //
    // Dump collection of points
    //
    friend std::ostream& operator << (std::ostream& os, PointsSpace & ps){

      PointId i = 0;
      BOOST_FOREACH(Points::value_type p, ps.points__){     
	os << "point["<<i++<<"]=" << p << std::endl;
      }
      return os;
    };

  public:

    PointsSpace(PointId num_points, Dimensions num_dimensions) 
      : num_points__(num_points), num_dimensions__(num_dimensions)
    {init_points();};

    inline const PointId getNumPoints() const {return num_points__;}
    inline const PointId getNumDimensions() const {return num_dimensions__;}
    inline const Point& getPoint(PointId pid) const { return points__[pid];}
 
  private:
    //
    // Init collection of points
    //
    void init_points();
    
    PointId num_points__;
    Dimensions num_dimensions__;
    Points points__;
  };

  // 
  //  This class represents a cluster
  // 
  class Clusters {

  private:
   
    ClusterId num_clusters__;    // number of clusters
    PointsSpace& ps__;           // the point space
    Dimensions num_dimensions__; // the dimensions of vectors
    PointId num_points__;        // total number of points
    ClustersToPoints clusters_to_points__;
    PointsToClusters points_to_clusters__;
    Centroids centroids__;

    //
    // Zero centroids
    //
    void zero_centroids();     

    //
    // Zero centroids
    //
    void compute_centroids();     
    
    //
    // Initial partition points among available clusters
    //
    void initial_partition_points();

  public:
   
    //
    // Dump ClustersToPoints
    //
    friend std::ostream& operator << (std::ostream& os, Clusters & cl){
      
      ClusterId cid = 0;
      BOOST_FOREACH(ClustersToPoints::value_type set, cl.clusters_to_points__){
	os << "Cluster["<<cid<<"]=(";
	BOOST_FOREACH(SetPoints::value_type pid, set){
	  Point p = cl.ps__.getPoint(pid);
	  os << "(" << p << ")";
	}
	os << ")" << std::endl;
	cid++;
      }
      return os;
    }
    

    Clusters(ClusterId num_clusters, PointsSpace & ps) 
      : num_clusters__(num_clusters), ps__(ps), 
	num_dimensions__(ps.getNumDimensions()),
	num_points__(ps.getNumPoints()),
	points_to_clusters__(num_points__, 0){

      ClusterId i = 0;
      Dimensions dim;
      for (; i < num_clusters; i++){
	Point point;   // each centroid is a point
	for (dim=0; dim<num_dimensions__; dim++) 
	  point.push_back(0.0);
	SetPoints set_of_points;

	// init centroids
	centroids__.push_back(point);  

	// init clusterId -> set of points
	clusters_to_points__.push_back(set_of_points);
	// init point <- cluster
      }
      /*
      std::cout << "Centroids" 
		<<std::endl<< centroids__;
      std::cout << "PointsToClusters" 
		<<std::endl<< points_to_clusters__;
      std::cout << "ClustersToPoints" 
		<<std::endl<< clusters_to_points__;
      */
    };
    
    //
    // k-means
    //
    void k_means (void);
  };

}





namespace Clustering{

  //
  // Dump a point
  //
  std::ostream& operator << (std::ostream& os, Point& p){
    
    BOOST_FOREACH(Point::value_type d, p){ os << d << " "; }    
    return os;
  }

  //
  // distance between two points
  //
  Distance distance(const Point& x, const Point& y)
  {
    Distance total = 0.0;
    Distance diff;
    
    Point::const_iterator cpx=x.begin(); 
    Point::const_iterator cpy=y.begin();
    Point::const_iterator cpx_end=x.end();
    for(;cpx!=cpx_end;++cpx,++cpy){
      diff = *cpx - *cpy;
      total += (diff * diff); 
    }
    return total;  // no need to take sqrt, which is monotonic
  }

  //
  // Dump collection of Points
  //
  std::ostream& operator << (std::ostream& os, Points& cps){
    BOOST_FOREACH(Points::value_type p, cps){ os<<p <<std::endl;}
    return os;
  }

  //
  // Dump a Set of points
  //
  std::ostream& operator << (std::ostream& os, SetPoints & sp){
    
    BOOST_FOREACH(SetPoints::value_type pid, sp){     
      os << "pid=" << pid << " " ;
    }
    return os;
  }

  //
  // Dump ClustersToPoints
  //
  std::ostream& operator << (std::ostream& os, ClustersToPoints & cp){
    ClusterId cid = 0;
    BOOST_FOREACH(ClustersToPoints::value_type set, cp){
      os << "clusterid["  << cid << "]" << "=(" 
	 << set << ")" << std::endl; 
      cid++;
    }
    return os;
  }

  //
  // Dump ClustersToPoints
  //
  std::ostream& operator << (std::ostream& os, PointsToClusters & pc){
    PointId pid = 0;
    BOOST_FOREACH(PointsToClusters::value_type cid, pc){
      
      std::cout << "pid[" << pid << "]=" << cid << std::endl;      
      pid ++;
    }
    return os;
  }

  //
  // Init collection of points
  //  
  void PointsSpace::init_points() { 
  			int i = 0;
		
		for (row_iterator it_a = (*a).begin(); it_a != (*a).end(); ++it_a, ++i)
			{
				
				doc>>idf; 
				*it_a =idf;
				//cout<<idf<<" ";
				
			}
			cout<<"\n A matrix number time iteration is: "<<i<<endl;
			
		i=0;
		//k-nn matrix initialization
		for(row_iterator2 it_b=(*b).begin();it_b!=(*b).end();++it_b,++i)
			{
				if(i%(k+1)==0)
					*it_b=i/(k+1);
				else
					*it_b=99999999;
			}
		
	    
	    cout<<"\n After k-nn matrix initialization";
		int j,p,d;
		double tempV1,tempV2;
		double vect1[dimension],vect2[dimension];
		row_iterator it_c=(*a).begin();
		row_iterator2 it_d=(*b).begin();
		
		double *distanceVect=(double*)malloc(k*sizeof(double));
		double *neighbors=(double*)malloc(k*sizeof(double));
		for(i=0;i < no_of_vector;i++)
			{
				++it_d;
				for(int ppp=0;ppp<k;ppp++)
				{
					distanceVect[ppp]=999999;
					neighbors[ppp]=99999999;
				}
				for(j=0;j<dimension;j++)
					{
						vect1[j]=*it_c;
						++it_c;
					}
				row_iterator it_a=(*a).begin();
					
				for(p=0;p<no_of_vector;p++)
					{
						
						for(j=0;j<dimension;j++)
							{
								vect2[j]=*it_a;
								++it_a;
							}
						/*double dist=distance2(vect1,vect2,dimension);
					
						if(distanceVect[k-1]>dist)
							{	
								distanceVect[k-1]=dist;
								neighbors[k-1]=p;
								d=k-1;
								while(d>0 && distanceVect[d-1]>distanceVect[d])
									{
										tempV2=distanceVect[d-1];
										tempV1=neighbors[d-1];
										distanceVect[d-1]=distanceVect[d];
										distanceVect[d]=tempV2;
										neighbors[d-1]=neighbors[d];
										neighbors[d]=tempV1;
										d--;	
									}	
					
							}*/
			
						}
						
				/*for(p=0;p<k;p++)
					{
						*it_d=neighbors[p];
						++it_d;
					}
				*/
		
			}
		
   
    for (PointId i=0; i < num_points__; i++){
      Point p;
      for (Dimensions d=0 ; d < num_dimensions__; d++)
	{ 
		
		p.push_back( vect1[d] );
	 }     
      points__.push_back(p);

      
    }
  }  

  //
  // Zero centroids
  //
  void Clusters::zero_centroids() {
      
    BOOST_FOREACH(Centroids::value_type& centroid, centroids__){
      BOOST_FOREACH(Point::value_type& d, centroid){
	d = 0.0;
      }
    }
  }

  //
  // Compute Centroids
  //
  void Clusters::compute_centroids() {

    
    Dimensions i;
    ClusterId cid = 0;
    PointId num_points_in_cluster;
    // For each centroid
    BOOST_FOREACH(Centroids::value_type& centroid, centroids__){

      num_points_in_cluster = 0;

      // For earch PointId in this set
      BOOST_FOREACH(SetPoints::value_type pid, 
		    clusters_to_points__[cid]){
	
	Point p = ps__.getPoint(pid);
	//std::cout << "(" << p << ")";
	for (i=0; i<num_dimensions__; i++)
	  centroid[i] += p[i];	
	num_points_in_cluster++;
      }
      //
      // if no point in the clusters, this goes to inf (correct!)
      //
      for (i=0; i<num_dimensions__; i++)
	centroid[i] /= num_points_in_cluster;	
      cid++;
    }
  }

  //
  // Initial partition points among available clusters
  //
  void Clusters::initial_partition_points(){
    
    ClusterId cid;
    
    for (PointId pid = 0; pid < ps__.getNumPoints(); pid++){
      
      cid = pid % num_clusters__;

      points_to_clusters__[pid] = cid;
      clusters_to_points__[cid].insert(pid);
    }
    
    /*
    std::cout << "Points_to_clusters " << std::endl;
    std::cout << points_to_clusters__;
    std::cout << "Clusters_to_points " << std::endl;
    std::cout << clusters_to_points__;
    */

  }

  //
  // k-means
  //
  void Clusters::k_means(void){
    
    bool move;
    bool some_point_is_moving = true;
    unsigned int num_iterations = 0;
    PointId pid;
    ClusterId cid, to_cluster;
    Distance d, min;
    

    //
    // Initial partition of points
    //
    initial_partition_points();

    //
    // Until not converge
    //
    while (some_point_is_moving){

      //std::cout << std::endl << "*** Num Iterations " << num_iterations  << std::endl << std::endl ;;

      some_point_is_moving = false;

      compute_centroids();
      //      std::cout << "Centroids" << std::endl << centroids__;      

      //
      // for each point
      //
      for (pid=0; pid<num_points__; pid++){
	
	// distance from current cluster
	min = 
	  Clustering::distance(centroids__[points_to_clusters__[pid]], 
			       ps__.getPoint(pid));

 	//std::cout << "pid[" << pid << "] in cluster=" << points_to_clusters__[pid] << " with distance=" << min << std::endl;

	//
	// foreach centroid
	//
	cid = 0; 
	move = false;
	BOOST_FOREACH(Centroids::value_type c, centroids__){
	  

	  d = Clustering::distance(c, ps__.getPoint(pid));
	  if (d < min){
	    min = d;
	    move = true;
	    to_cluster = cid;

	    // remove from current cluster
	    clusters_to_points__[points_to_clusters__[pid]].erase(pid);

	    some_point_is_moving = true;
	   // std::cout << "\tcluster=" << cid << " closer, dist=" << d << std::endl;	    
	  }
	  cid++;
	}
	
	//
	// move towards a closer centroid 
	//
	if (move){
	  
	  // insert
	  points_to_clusters__[pid] = to_cluster;
	  clusters_to_points__[to_cluster].insert(pid);
	  //std::cout << "\t\tmove to cluster=" << to_cluster << std::endl;
	}
      }      

      num_iterations++;
    } // end while (some_point_is_moving)

    //std::cout << std::endl << "Final clusters" << std::endl;
    //std::cout << clusters_to_points__;
  }

  
}
using namespace Clustering;

int main()
{
  
  ClusterId num_clusters =7; 
  PointId num_points = size;
  Dimensions num_dimensions =dimension;
  doc.open("/media/livingroom/stxxl/stxxl-1.4.1/local/myfile_final64.txt");
  matrix_operation_statistic_data matrix_stats_before, matrix_stats_after;
	    matrix_stats_before.set();
	    stats_data stats_before = *stats::get_instance();
  stxxl::stats_data begin_1(*stats::get_instance());	
			
  PointsSpace ps(num_points, num_dimensions);
  

  Clusters clusters(num_clusters, ps);

  clusters.k_means();
  	stxxl::stats_data begin_2(*stats::get_instance());
	          	
	STXXL_MSG(begin_2-begin_1);
	
	stats_data stats_after = *stats::get_instance();
	matrix_stats_after.set();
	STXXL_MSG(stats_after - stats_before);
		
  
  doc.close();
  return 0;

}
