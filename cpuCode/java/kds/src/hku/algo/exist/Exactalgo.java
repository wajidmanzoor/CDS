package hku.algo.exist;

import java.util.Arrays;
import java.util.Map;
import java.util.Map.Entry;
import hku.algo.maxflow.FindMinCut;

public class Exactalgo {
	
	/** data structure used to save all matchings */
	private Map<String,int[]> Motif_Record=null;
	/** the number of vertex in the given motif*/
	private int motif_size=0;
	/** the number of vertex in the given graph*/
	private int graph_size=0;
	/** vertex's motif-degree*/
	private int[] Motif_degree=null;
	
	/**
	 * 
	 * @param map data structure used to record all matchings.
	 * @param motif_size the number of vertices in the given motif
	 * @param graph_size the number of vertices in the given graph
	 * @param Motif_degree the number of participating motifs per vertex
	 */
	public Exactalgo(Map<String,int[]> map,int motif_size,
			int graph_size,int[] Motif_degree) {
		this.Motif_Record=map;
		this.motif_size=motif_size;
		this.graph_size=graph_size;
		this.Motif_degree=Motif_degree;
	}
	
	/**
	 * 
	 * @param l
	 * @param u
	 * @param motif_num
	 * @return the min (S-T) cut
	 */
	public int[] Exact(double l,double u,long motif_num) {
		
		//wm: construct a flow network 
		FlowNetwork flownetwork=new FlowNetwork(Motif_Record, motif_size, graph_size, Motif_degree);

		//wm: guessed motif density 
		double alph=(u+l)/2;

		//wm: minimum motif density difference
		double bais=1.0/(graph_size*(graph_size-1));

		//wm: stopping condition 
		if(bais<0.000000000000001){
			bais=0.000000000000001;
		}

		//wm: generate a flownetwrok network 
		// return the constructed flow network as map of hash maps (index start of edge, key end of edge, values array (capacity, remaining flow ))
		Map<Integer,double[]>[] Network=flownetwork.Construct(alph);
		

		//wm: stores the S-T min cut
		FindMinCut compute=new FindMinCut(Network, Network.length-2, Network.length-1);

		//wm: store the result flow and average density
		double res_flow=0,res_alph=0;

		//wm : 
		int res[]=new int[graph_size];
		Arrays.fill(res, 1);
		

		
		while(u-l>bais) {
			double temp=compute.EdmondsKarp();

			
			if(temp==motif_num*motif_size) {

				//wm: when S == {s}. set the upper limit to the guess
				// which means subgraph has at max guess motif density 
				u=alph;
				

			}else {

				//wm: id S != {s}, verticies other than s are in S set of S-T cut
				//wm: means the subgraph S-{s} has atleast alph motif density 

				//wm: set lower bound to guess
				l=alph;


				res_alph=alph;
				res_flow=temp;

				//wm: get the densest subgraph 
				int temp_array[]=compute.getparent();
				for(int i=0;i<graph_size;++i) {
					res[i]=temp_array[i];
				}
				

			}

			//wm: set new guess
			alph=(u+l)/2;

			//wm: update the network, by updating the capacity and availble flow of edges 
			Network=flownetwork.Update(alph);
			
		}
		
		//wm: return the densest subgraph of the connected component
		return res;
	}
	
	
public int[] Exact(double l,double u,long motif_num,long n2) {
		
		FlowNetwork flownetwork=new FlowNetwork(Motif_Record, motif_size, graph_size, Motif_degree);
		double alph=(u+l)/2;
		double bais=1.0/(n2*(n2-1));
		if(bais<0.000000000000001){
			bais=0.000000000000001;
		}
		Map<Integer,double[]>[] Network=flownetwork.Construct(alph);
		
		FindMinCut compute=new FindMinCut(Network, Network.length-2, Network.length-1);
		//System.out.println("eeee");
		double res_flow=0,res_alph=0;
		int res[]=new int[graph_size];
		Arrays.fill(res, 1);
		
//		int a=0;
//		double temp_r=0;
		
		while(u-l>bais) {
			double temp=compute.EdmondsKarp();
//			temp=temp+temp_r;
			System.out.println(u+" "+l+" "+alph+" "+temp+" "+motif_num*motif_size);

			
			if(temp==motif_num*motif_size) {
				u=alph;
				
//				alph=(u+l)/2;				
//				Network=flownetwork.Update(alph);
//				compute.FlowNetwork=Network;
//				temp_r=0;
				
			}else {
//				if(a==0){
//					flownetwork.copyNetwork();
//				}
//				a++;
				l=alph;
				res_alph=alph;
				res_flow=temp;
				int temp_array[]=compute.getparent();
				for(int i=0;i<graph_size;++i) {
					res[i]=temp_array[i];
				}
				
			}
			alph=(u+l)/2;
			Network=flownetwork.Update(alph);
			
		}
		

		return res;
	}
	
	public boolean Try(double l,long motif_num) {
		FlowNetwork flownetwork=new FlowNetwork(Motif_Record, motif_size, graph_size, Motif_degree);
		double alph=l;
		//double bais=1.0/(graph_size*(graph_size-1));
		Map<Integer,double[]>[] Network=flownetwork.Construct(alph);
		FindMinCut compute=new FindMinCut(Network, Network.length-2, Network.length-1);
		//double res_flow=0,res_alph=0;
		//int res[]=new int[graph_size];
		double temp=compute.EdmondsKarp();
		//System.out.println(u+" "+l+" "+alph);
		if(temp==motif_num*motif_size&&temp!=graph_size*alph*motif_size) {
			return false;
		}else {
			return true;
		}
	}

}
