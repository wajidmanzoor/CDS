package hku.algo.prune;

import java.util.Arrays;

public class LocateCore {
	
	//wm: graph adjacency list
	private int[][] Graph=null;

	//wm: verticies, motif degree, average motif density, core values 
	private double[][] core=null;
	private int graph_size;
	public int[] reverse_map = null;
	
	public LocateCore(int[][] Graph,double[][] core,int graph_size) {
		this.Graph=Graph;
		this.core=core;
		this.graph_size=graph_size;
	}
	
	public DensestCore locate() {
		//wm: find the low bound, based on maximum average motif density, take the lower bound as the motif degree of that vertex
		double max=core[0][2];
		double kmax=0;
		for(int i=1;i<graph_size;++i) {
			if(max<core[i][2])
				max=core[i][2];
			if(kmax<core[i][1])
				kmax=core[i][1];
		}
			
		int low_bound=(int) Math.ceil(max);

		
		//wm: delete all verticies with motif degree lower than the lower bound 
		// PROBLEM: WE SHOULD CHECK CORE[INDEX][2]
		int index=1;
		int delete[]=new int[graph_size];
		Arrays.fill(delete, 0);
		for(;index<graph_size;++index) {
			if(core[index][1]>=low_bound)
				break;

			//wm: set verticies that needs to be deleted to -1
			delete[(int) core[index][0]]=-1;
		}

		int temp=0;
		for(int i=0;i<graph_size;++i) {
			if(delete[i]==0) {

				//wm: vertices that are not deleted
				delete[i]=temp;
				temp++;
			}
		}
		
		//wm: stores new vertices, motif degree, average motif density 
		double[][] new_core=new double[temp][2];
		

		//wm: store new graph
		int New_graph_size=temp;
		int New_Graph[][]=new int[New_graph_size][];
		temp=0;

		//wm: update the vertices in new core
		for(int i=index;i<graph_size;++i) {
			int m=(int)core[i][0];
			new_core[temp][0]=delete[m];
			new_core[temp][1]=core[i][1];
			temp++;
		}
		
		
		
		//wm: strore new graph as adjacency list
		for(int i=0;i<graph_size;++i) {

			//wm: vertices that are not deleted
			if(delete[i]!=-1) {
				int count=0;
				int array[]=new int[Graph[i].length];

				//wm: get neighbors that are not deleted and their count
				for(int j=0;j<Graph[i].length;++j) {
					if(delete[Graph[i][j]]!=-1) {
						array[count]=delete[Graph[i][j]];
						count++;
					}					
				}

				//wm: add neighbors that are not deleted
				New_Graph[delete[i]]=new int[count];
				for(int j=0;j<count;++j) {
					New_Graph[delete[i]][j]=array[j];
				}
				
			}
		}

		reverse_map = new int[graph_size];  // size is the number of surviving vertices

		for (int i = 0; i < delete.length; i++) {
			if (delete[i] != -1) {
				reverse_map[delete[i]] = i;
			}
		}


		//wm: number of removed motifs core[0][3] : total motifs, core[index-1][3]: total remaining motifs after deletion
		int delete_motif=(int)(core[0][3]-core[index-1][3]);

		int motif_degree[]=new int[New_graph_size];
		
		DensestCore result=new DensestCore(New_Graph, New_graph_size, low_bound, index-1, delete_motif,core[index-1][2],(int)kmax);
		return result;
	}
	

}
