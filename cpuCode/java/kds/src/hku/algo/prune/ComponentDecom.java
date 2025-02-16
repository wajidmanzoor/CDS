package hku.algo.prune;

import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;

public class ComponentDecom {
	
	private int[][] Graph=null;
	private int graph_size;
	private Map<String,int[]> motif_list;
	
	public ComponentDecom(int[][] Graph,int graph_size,Map<String,int[]> motif_list) {
		this.Graph=Graph;
		this.graph_size=graph_size;
		this.motif_list=motif_list;
	}
	
	public Queue<Component> decompose() {
		int delete[]=new int[graph_size];
		Arrays.fill(delete, 0);
		
		//wm: perform a BFS to find all connected components, vertices in each component has the same value in delete array. 
		int index=0;
		for(int i=0;i<graph_size;++i) {
			if(delete[i]==0) {
				index++;
				BFS(delete,i,index);				
			}
		}

		//wm : queue of connected components 
		Queue<Component> result=new LinkedList<Component>();
		
		//wm: if only one connected conponents i.e whole graph is connected
		if(index==1) {
			long motif_num=0;
			int motif_degree[]=new int[graph_size];
			int motif_size=1;
			for(Entry<String,int[]> entry:motif_list.entrySet()) {
				int temp[]=entry.getValue();
				for(int i=0;i<temp.length-1;++i) {
					motif_degree[temp[i]]+=temp[entry.getValue().length-1];
				}
				motif_size=entry.getValue().length-1;
				motif_num+=temp[motif_size];
				
			}
			motif_num=motif_num;
			Component c=new Component(Graph,graph_size,
					motif_list,motif_num,(double)motif_num/(graph_size*1.0),motif_degree);
			result.add(c);
		
		//wm: if more than one connected conponents.
		}else {
			

			//wm: array that stores maps of motifs in each component, index 0 stores full, index 1 stores motiffs of connected component 1 
			Map<String,int[]>[] map_array=new Map[index+1];
			for(int i=1;i<index+1;++i) {
				map_array[i]=new HashMap<String,int[]>();
			}
			
			
			//wm: stores graph size of connected components
			int New_graph_size[]=new int[index+1];

			//wm: stores mapping from old verticies to new vertices of connected components
			int Map_node[]=new int[graph_size];


			Arrays.fill(Map_node, 0);
			Arrays.fill(New_graph_size, 0);
			
			//wm: iter through each vertex to form the new mapping of verticies and graph size of each connected component
			for(int i=0;i<graph_size;++i) {
				Map_node[i]=New_graph_size[delete[i]];
				New_graph_size[delete[i]]++;
			}
			

			//wm: store the updated mottifs in new array of hash maps (as vertex numbering was changed)
			for(Entry<String,int[]> entry:motif_list.entrySet()) {
				int temp=entry.getValue()[0];
				int array[]=entry.getValue();
				for(int i=0;i<array.length-1;++i) {
					array[i]=Map_node[array[i]];
				}
				map_array[delete[temp]].put(entry.getKey(),entry.getValue());
			}
			
			
			//wm: store adjancey list for each new connected component
			//index 1 stores adjancey list of component 1 
			int [][][] C_Graph=new int[index+1][][];

			//wm: construct the graph vertex array for each connected component
			for(int i=1;i<index+1;++i) {
				C_Graph[i]=new int[New_graph_size[i]][];
			
			}

			//wm: get adjancey list of each vertex (new numbering) in each connected component
			for(int i=0;i<graph_size;++i) {

				//wm: stores neighbors 
				int temp[]=new int[Graph[i].length];
				for(int j=0;j<Graph[i].length;++j) {
					temp[j]=Map_node[Graph[i][j]];
				}

				//wm: stores neighbors of vertex i (map_node[i]) of connected component delete[i] 
				C_Graph[delete[i]][Map_node[i]]=temp;
			}


			//wm: update the motif degree for each component and add to  connected component queue
			int motif_size=1;
			for(int i=1;i<index+1;++i) {
				long motif_num=0;

				//wm: motif degree array
				int motif_degree[]=new int[New_graph_size[i]];
				Arrays.fill(motif_degree, 0);

				//wm: iter through all motifs 
				for(Entry<String,int[]> entry:map_array[i].entrySet()) {

					//wm: increament motif count 
					motif_num+=entry.getValue()[entry.getValue().length-1];
					motif_size=entry.getValue().length-1;

					//wm: store vertices in the motif
					int temp[]=entry.getValue();
					for(int j=0;j<temp.length-1;++j)

						//wm: update the degree of vertex by 1
						motif_degree[temp[j]]+=temp[temp.length-1];
				}
				
				//wm: create connected component object 
				Component c=new Component(C_Graph[i],New_graph_size[i],
						map_array[i],motif_num,(double)motif_num/(New_graph_size[i]*1.0),motif_degree);
				
				//wm: add to queue
				result.add(c);
			}
		}
		return result;
	}
	
	private void BFS(int[] delete,int s,int index) {
		Queue<Integer> queue=new LinkedList<Integer>();
		queue.add(s);
		while(!queue.isEmpty()) {
			int node=queue.poll();
			delete[node]=index;
			for(int i=0;i<Graph[node].length;++i) {
				if(delete[Graph[node][i]]==0) {
					queue.add(Graph[node][i]);
				}
			}
		}
	}
}
