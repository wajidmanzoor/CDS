package hku.algo.prune;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

public class InvalidEdgePruning {
	
	private Map<String,int[]> motif_list=null;
	
	private int[][] Graph=null;
	
	private int graph_size;
	
	public InvalidEdgePruning(Map<String,int[]> motif_list,
			int[][] Graph,int graph_size) {

		//wm: map of motif : vertices of motif
		this.motif_list=motif_list;
		this.Graph=Graph;
		this.graph_size=graph_size;
	}
	
	public int Prune() {


		int count=0;

		//wm: list of hash maps
		Map<Integer,Integer> neig[]=new Map[graph_size];

		//wm: each hashmap stores the edges that are part of the motif, that is add valid neighbors
		for(int i=0;i<graph_size;++i) {
			neig[i]=new HashMap<Integer,Integer>();
		}

		//wm: iter over all the motifs 
		for(Entry<String,int[]> entry:motif_list.entrySet()) {
			int temp[]=entry.getValue();

			//wm: iter over the vertices in the motif
			for(int i=0;i<temp.length-1;++i) {
				for(int j=0;j<temp.length-1;++j) {

					//wm: get all edges except circular ones 
					if(i!=j) {
						if(!neig[temp[i]].containsKey(temp[j])) {

							//wm: add valid neighbors in the neighbor list
							neig[temp[i]].put(temp[j], 0);
						}
					}
				}
			}
		}


		//wm: creae new graph that only include the edges that are part of motiffs 
		for(int i=0;i<graph_size;++i) {
			int temp[]=new int[Graph[i].length];
			
			//wm: get the new size of neighbors
			int size=0;


			for(int j=0;j<Graph[i].length;++j) {
				if(neig[i].containsKey(Graph[i][j])) {
					temp[size]=Graph[i][j];
					size++;
				}
			}

			//wm: get new neighbor list 
			if(size<Graph[i].length) {
				int array[]=new int[size];
				for(int j=0;j<size;++j)
					array[j]=temp[j];

				//wm: number of removed edges
				count+=(Graph[i].length-size);
				Graph[i]=array;
				
			}
		}

		//wm: divided by as each edge will be considered twice
		return count/2;
	}

}
