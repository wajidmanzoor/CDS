package hku.algo.cds;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;

import hku.algo.findgeneralpattern.FindMotif;
import hku.util.Log;

/**
 * 
 * @author yukaiqiang
 *
 */
public class CDSdecompose {

	/** adjacent list of the given graph */
	private int[][] Graph = null;
	/** adjacent matrix of the given motif */
	private int[][] Motif = null;
	/** the number of vertex in the given graph */
	private int graph_size;
	/** the number of vertex in the given motif */
	private int motif_size;
	/** the degree for each vertex in the motif */
	private long[] degree;
	/** data structure to record the information of sharing motifs */
	private Map<Integer, Integer> Share[] = null;
	private int motif_type=1;
	
	private Map<String,int[]> motif_list=null;
	
	private long[] motif_degree;

	/**
	 * 
	 * @param Graph
	 *            adjacent list of the given graph
	 * @param Motif
	 *            adjacent matrix of the given motif
	 * @param graph_size
	 *            the number of vertex in the given graph
	 * @param motif_size
	 *            the number of vertex in the given motif
	 */
	public CDSdecompose(int[][] Graph, int[][] Motif, int graph_size, int motif_size,int motif_type,
			Map<String,int[]> motif_list,long[] motif_degree) {
		this.motif_list=motif_list;

		//wm: adjacent list
		this.Graph = Graph;
		//wm: adjacent matrix
		this.Motif = Motif;

		this.graph_size = graph_size;
		this.motif_size = motif_size;
		this.motif_type=motif_type;
		this.motif_degree=motif_degree;

		//wm: get degree of each vertex in motif from the adjancey matrix
		degree = new long[motif_size];
		for (int i = 0; i < motif_size; ++i) {
			int count = 0;
			for (int j = 0; j < motif_size; ++j) {
				if (Motif[i][j] > 0)
					count++;
			}
			degree[i] = count;
		}

		Share = new Map[graph_size];
		for (int i = 0; i < graph_size; ++i) {
			Share[i] = new HashMap<Integer, Integer>();
		}

	}


	/**
	 * 
	 * @return
	 */
	public double[][] Decompose() {
		
		long max = 0;
		int motif_num = 0;
		int count = 0;

		int mark[]=new int[graph_size];
		Arrays.fill(mark, 0);
		Map<Integer, Long> map;
		
		int array_index[]=new int[graph_size];
		Arrays.fill(array_index, 0);
		long array_degree[]=new long[graph_size];

		int new_array[]=new int[graph_size];
		Arrays.fill(new_array, 0);
		int new_map[]=new int[graph_size];
		Arrays.fill(new_map, 0);

		//wm: Graph and motif arrays
		KList k=new KList(Graph,motif_size);

		//wm: Lists all motif_size- cliques of the graph, to get the motif degree 
		k.ListFast();

		//wm:  Gives the number of motiff a vertex is part of
		array_degree=k.getMotifDegree();

		motif_degree=array_degree;

		System.out.println("Motif Degree");
		for(int i=0;i<graph_size;i++){
			System.out.print(motif_degree[i]+" ");
		}
		System.out.println();

		//wm: copy of motif degree
		long temp_degree[]=new long[graph_size];
		for(int i=0;i<graph_size;++i) {
			temp_degree[i]=motif_degree[i];
		}
		
		//wm: get max degree and average motif degree degree 
		for (int i = 0; i < graph_size; ++i) {
			if (max < motif_degree[i])
				max = motif_degree[i];
			motif_num += motif_degree[i];
		}


		max=max+1;

		//wm: gives actual motif number as each motif might have been calculated motif size time ???
		motif_num = motif_num / motif_size;

		// data structure used to save the result
		double[][] result = new double[graph_size+1][5];
		result[0][2] = motif_num / (double) graph_size;
		result[0][3] =motif_num;


	
		
		
		/****
		 * debug 
		 */
		if(false){
			//long bin[]=new long[max+1];
			int multi=(int) ((max+1)/Integer.MAX_VALUE+1);
			int bin[][]=new int[multi][];
			for(int i=0;i<multi;++i){
				bin[i]=new int[Integer.MAX_VALUE];
			}
			for(int i=0;i<multi;++i)
				Arrays.fill(bin[i], 0);
			for(int i=1;i<graph_size+1;++i) {
				int a=(int) (motif_degree[i-1]/Integer.MAX_VALUE);
				int b=(int) (motif_degree[i-1]%Integer.MAX_VALUE);
				bin[a][b]+=1;
			}
			
			int start=1;
			for(long d=0;d<=max;++d) {
				int a=(int) (d/Integer.MAX_VALUE);
				int b=(int) (d%Integer.MAX_VALUE);
				int num=bin[a][b];
				bin[a][b]=start;
				start+=num;
			}
			
			int pos[]=new int[graph_size+1];
			int vert[]=new int[graph_size+1];
			for(int v=1;v<graph_size+1;++v) {
				int a=(int) (motif_degree[v-1]/Integer.MAX_VALUE);
				int b=(int) (motif_degree[v-1]%Integer.MAX_VALUE);
				pos[v]=bin[a][b];
				vert[pos[v]]=v;
				bin[a][b]+=1;
			}
			
			for(long d=max;d>=1;d--) {
				int a=(int) (d/Integer.MAX_VALUE);
				int b=(int) (d%Integer.MAX_VALUE);
				int c=(int) ((d-1)/Integer.MAX_VALUE);
				int d1=(int) ((d-1)%Integer.MAX_VALUE);
				bin[a][b]=bin[c][d1];
			}
			bin[0][0]=1;
			long core_number=0;
			
			for (int i = 1; i < graph_size+1; ++i) {
				
				
				int index=vert[i]-1;
				

				// get the vertex with minimum motif-degree
//				if (update_min < index_min) {
//					index_min = update_min;
//				} else {
//					for (; index_min < max; ++index_min) {
//						if (!node_list[index_min].isEmpty())
//							break;
//					}
//				}
//				update_min = Integer.MAX_VALUE;
//				int index = (Integer) node_list[index_min].remove(0);

				// record the result
				
				count++;
				result[count][0] = index;
				result[count][1] = motif_degree[index];
				motif_num -= temp_degree[index];
				if(motif_num<0) motif_num=0;
				result[count][3]=motif_num;
//				
				
//				if(core_number>motif_degree[index]) {
//					System.out.println(count+" break");
//					return result;
//				}motif_num -= temp_degree[index];
				
				if (graph_size - count > 0) {
					result[count][2] = motif_num / (double) (graph_size - count);
				} else {
					result[count][2] = 0;
				}
//					
				core_number=motif_degree[index];
//				if(core_number>motif_degree[index]) {
//					result[count][4]=core_number;
//				}else {
//					core_number=motif_degree
//				}
				

				// update
				if(motif_degree[index]>0) {
					
					
					//map=Update(index,mark);
//					System.out.println("****0"+" "+index);
//					for(Entry<Integer, Integer> entry : map.entrySet()) {
//						System.out.println(entry.getKey()+" "+entry.getValue());
//					}
//					System.out.println("*****1"+" "+mark[37]);
					//System.out.println(111);
					map=Generate(index,mark,new_array,new_map);
					//System.out.println(222);
//					for(Entry<Integer, Integer> entry : map.entrySet()) {
//						System.out.println(entry.getKey()+" "+entry.getValue());
//					}
//					System.out.println("*****2");
					if(!map.isEmpty())
					for(Entry<Integer, Long> entry : map.entrySet()) {
						int temp_key=entry.getKey();
						Long temp_value=entry.getValue();
						//node_list[motif_degree[temp_key]].remove((Integer)temp_key);			
						if(motif_degree[temp_key]>motif_degree[index]) {
							long du=motif_degree[temp_key];
							int pu=pos[temp_key+1];
							int a=(int) (du/Integer.MAX_VALUE);
							int b=(int) (du%Integer.MAX_VALUE);
							int pw=bin[a][b];
							int w=vert[pw];
							if((temp_key+1)!=w) {
								pos[temp_key+1]=pw;
								vert[pu]=w;
								pos[w]=pu;
								vert[pw]=temp_key+1;							
							}
							bin[a][b]+=1;
							
						}
						motif_degree[temp_key]-=temp_value;
						temp_degree[temp_key]-=temp_value;
						
//						if(update_min>motif_degree[temp_key]) {
//							update_min=motif_degree[temp_key];
//						}
//						node_list[motif_degree[temp_key]].add(temp_key);
					}
				}else {
					motif_degree[index]=0;
					//System.out.println(".");
				}
					
				//node_list[index_min].remove((Integer)index);
				mark[index]=1;
				
				

			}
		}else{
			
		
			// wm: get count of vertcies who has a particular motif degree, eg motif degree 0 - 3 verticies.
			long bin[]=new long[(int) (max+1)];
			Arrays.fill(bin, 0);
			for(int i=0;i<graph_size;++i) {
				bin[(int) motif_degree[i]]+=1;
			}
			
			//wm: cumulative prefix sum 
			Long start=0L;
			for(int d=0;d<=max;++d) {
				long num=bin[d];
				bin[d]=start;
				start+=num;
			}
			
		    
			//wm: sorts the vertex by degree asscending order (kinda like bucket sort)
			
     		//wm: stores position of verticies in sorted array.
			int pos[]=new int[graph_size+1];

		    //wm: stores sorted array of verticies by motif degree 
			int vert[]=new int[graph_size+1];
			for(int v=0;v<graph_size;++v) {
				pos[v]=(int) bin[(int) motif_degree[v]];
				vert[pos[v]]=v;
				bin[(int) motif_degree[v]]+=1;
			}


		    //wm: reset the bin to before sort state by shifting one element to right
			for(long d=max;d>=1;d--) {
				bin[(int) d]=bin[(int) (d-1)];
			}

			//idk: wm: dont know why ?
			if(bin.length!=0)
				bin[0]=1;
			long core_number=0;
			//long aaaaa=0;
			
			for (int i = 0; i < graph_size; ++i) {
				
				//wm: position of vertex with min motif degree
				int index=0;

				//wm: min motif degree
				long index_min=0XFFFFFF;
				//wm: get the vertex with smallest motif degree which isn't processed yet 
				for(int yuu=0;yuu<motif_degree.length;++yuu) {
					if(index_min>motif_degree[yuu]&&mark[yuu]==0) {
						index=yuu;
						index_min=motif_degree[yuu];
					}
				}
				
				System.out.println("Removed " +index);

				
				
				count++;
				result[count][0] = index;
				result[count][1] = motif_degree[index];
				if(motif_degree[index]>0) {

					//wm: extract a neighborhood arround the vertex which is just 2 edges away and get the motif degree of that neighborhood
					//wm: here the motif degree is only cacluated using the motifs that include the desired vertex (index)
					//map vertex  of subgrah : motif degree 
					map=Generate(index,mark,new_array,new_map);

					//wm: motif degree that will be excluded if the vertex is removed
					long delete_count=0;
					if(!map.isEmpty())
					for(Entry<Integer, Long> entry : map.entrySet()) {

						//wm: vertex 
						int temp_key=entry.getKey();

						//wm: motif degree
						Long temp_value=entry.getValue();

						//wm: motif degree that will be excluded if the vertex is removed
						delete_count+=temp_value;

						//wm: decrease motif degree of the vertex 
						motif_degree[temp_key]-=temp_value;
						temp_degree[temp_key]-=temp_value;
					}
					delete_count=delete_count/(motif_size-1);
					
					motif_num -= delete_count;
					
					result[count][3]=motif_num;
					if (graph_size - count > 0) {
						result[count][2] = motif_num / (double) (graph_size - count);
					} else {
						result[count][2] = 0;
					}
					
				}else {
					result[count][3]=motif_num;
					if (graph_size - count > 0) {
						result[count][2] = motif_num / (double) (graph_size - count);
					} else {
						result[count][2] = 0;
					}
					motif_degree[index]=0;
				}
				
				
				//wm: mark vertex as processed
				mark[index]=1;
				
				

			}
		}
		
		/**
		 result[i][0]	Vertex removed in step i
		result[i][1]	Motif-degree of removed vertex before deletion
		result[i][2]	motif density: Average motif count per remaining vertex (motif_num / (graph_size - count)) 
		result[i][3]	Remaining total motifs in the graph after removal
		result[i][4]	Core number (motif-core value) of the removed vertex
		 */
		return result;

	}
	
	public Map<Integer, Long> Generate(int index,int mark[],int array[],int map_s[]){
		LinkedList temp_list=new LinkedList<Integer>();

		//wm: add index to temp list 
		temp_list.add(index);

		array[index]=1;
		Queue queue=new LinkedList();
		queue.add(index);
		int d=1;

		//wm: do a DFS in neighborhood of vertex (which are <= 2 edges away)
		//wm: 
		while(!queue.isEmpty()&&d+1<=2) {
			int temp=(Integer) queue.poll();
			d=array[temp];
			for(int i=0;i<Graph[temp].length;++i) {

				//wm: if they are not processed and are only <=2 edges away , add them to temp_list
				if(array[Graph[temp][i]]==0&&d+1<=2&&mark[Graph[temp][i]]==0) {
					queue.add(Graph[temp][i]);
					array[Graph[temp][i]]=d+1;
					temp_list.add(Graph[temp][i]);
				}
			}
		}
		
		//wm: gives the number of verticies in neighboorhood (<= 2 edges away)
		int count=temp_list.size();


		int map_array[]=new int[count];
		int num=0;
		int new_graph[][]=new int[count][];


		//wm: create a adjancy list if the new subgraph 
		for(int i=0;i<count;++i) {
			int node=(Integer) temp_list.get(i);
			map_array[i]=node;
			map_s[node]=num;
			
				int temp_count=0;
				//wm: get count of neighbors that are processed
				for(int j=0;j<Graph[node].length;++j) {
					if(array[Graph[node][j]]!=0&&Graph[node][j]!=node) {
						temp_count++;
					}
				}
				new_graph[num]=new int[temp_count];
				temp_count=0;

				//wm: create adjancey list of new subgraph
				for(int j=0;j<Graph[node].length;++j) {
					if(array[Graph[node][j]]!=0&&Graph[node][j]!=node) {
						new_graph[num][temp_count]=Graph[node][j];
						temp_count++;
					}
				}
				num++;
		}

		//wm: update the status in array back to 0. 
		for(int i=0;i<count;++i) {
			int update=(Integer) temp_list.get(i);
			array[update]=0;
		}
			
		//wm: reindex the subgraph 
		for(int i=0;i<count;++i) {
			for(int j=0;j<new_graph[i].length;++j) {
				new_graph[i][j]=map_s[new_graph[i][j]];
			}
		}
		
		//wm : apply core decompose to get vertices sorted by core values (descending) for the subgraph 
        // also order, order of each vertex in sorted array of vertices by core values (descending)
		KList f=new KList(new_graph,motif_size);

		
		//wm: extract all k cliques that have the desired vertex and get the motif degree 
		f.ListOne(0);
		long[] t_a=f.getMotifDegree();
		Map<Integer,Long> result=new HashMap<Integer,Long>();
		for(int i=1;i<t_a.length;++i) {
			if(t_a[i]>0)
				result.put(map_array[i], t_a[i]);
				System.out.println(map_array[i]);
				System.out.println(t_a[i]);
				System.out.println();

		}


		
		
		// wm: motif degree of each vertex in subgraph that is because of a motif that includes the desired vertex
		return result;
	}
	
	
	public int com(int index,int mark[],int array[],int map_s[]){
		LinkedList temp_list=new LinkedList<Integer>();
		temp_list.add(index);
		//int[] array=new int[graph_size];
		//Arrays.fill(array, 0);
		array[index]=1;
		Queue queue=new LinkedList();
		queue.add(index);
		int d=1;
		while(!queue.isEmpty()&&d+1<=motif_size) {
			int temp=(Integer) queue.poll();
			d=array[temp];
			for(int i=0;i<Graph[temp].length;++i) {
				if(array[Graph[temp][i]]==0&&d+1<=motif_size) {
					queue.add(Graph[temp][i]);
					array[Graph[temp][i]]=d+1;
					temp_list.add(Graph[temp][i]);
				}
			}
		}
		
		int count=temp_list.size();
		int map_array[]=new int[count];	
		int num=0;
		int new_graph[][]=new int[count][];
		for(int i=0;i<count;++i) {
			int node=(Integer) temp_list.get(i);
			map_array[i]=node;
			map_s[node]=num;
			
				int temp_count=0;
				for(int j=0;j<Graph[node].length;++j) {
					if(array[Graph[node][j]]!=0&&Graph[node][j]!=node) {
						temp_count++;
					}
				}
				new_graph[num]=new int[temp_count];
				temp_count=0;
				for(int j=0;j<Graph[node].length;++j) {
					if(array[Graph[node][j]]!=0&&Graph[node][j]!=node) {
						new_graph[num][temp_count]=Graph[node][j];
						
						temp_count++;
					}
				}
				num++;
		}
		for(int i=0;i<count;++i) {
			int update=(Integer) temp_list.get(i);
			array[update]=0;
		}
		
//		System.out.println();
//		for(int i=0;i<count;++i) {
//			for(int j=0;j<new_graph[i].length;++j) {
//				System.out.print(new_graph[i][j]+" ");
//			}
//			System.out.println();
//		}
		//System.out.println("******");
		for(int i=0;i<count;++i) {
			//System.out.println();
			for(int j=0;j<new_graph[i].length;++j) {
				new_graph[i][j]=map_s[new_graph[i][j]];
				//System.out.print(new_graph[i][j]+" ");
			}
		}
		//System.out.println();
		FindMotif f=new FindMotif( Motif,new_graph,count,motif_type);
		int mmm=f.Match_V3();
		//System.out.println(mmm);
		//Log.write("111112");
		return mmm;
	}


}



//debug: replace DFS by BFS
//public int[] ComputeMD() {
//	// initialize the array
//	int motif_degree[] = new int[graph_size];
//	Arrays.fill(motif_degree, 0);
//
//	int H[] = new int[motif_size];
//	ArrayList record[] = new ArrayList[motif_size];
//	int D[] = new int[graph_size];
//	int F[] = new int[motif_size];
//	
//
//	for (int i = 0; i < graph_size; ++i) {
//		Arrays.fill(F, 0);
//		int motif_num = 0;
//		for (int j = 0; j < motif_size; ++j)
//			record[j] = new ArrayList<Integer>();
//		Arrays.fill(D, 0);
//		int deep = 0;
//		H[0] = i;
//		record[0].add(0);
//		while (deep >= 0) {
//			if (deep == motif_size - 1) {
//				deep--;
//				if(H[0]==0) {
//					System.out.println();
//					for(int mm=0;mm<motif_size;++mm)
//						System.out.print(H[mm]);
//					System.out.println();
//				}
//				int temp = CountMotif(H);
//				//System.out.println(temp);
//				motif_num += temp;
//				// Construct(H, temp);
//				continue;
//			}
//
//			int index = H[deep];
//			// mark that we have visited this vertex
//			
//			D[index] = 1;
//			// find a feasible vertex
//			int j = F[deep];
//			for (; j < Graph[index].length; ++j) {
//				if (D[Graph[index][j]] == 0)
//					break;
//			}
//			
//			if(H[0]==0) {
//				System.out.println("&&&&&");
//				for(int mm=0;mm<motif_size;++mm)
//				System.out.print(H[mm]+" ");
//				System.out.println();
//				System.out.println("(()())");
//			}
//			
//			if (j < Graph[index].length) {
//				F[deep] = j;
//				H[deep + 1] = Graph[index][j];
//				D[Graph[index][j]] = 1;
//				record[deep + 1].add(Graph[index][j]);
//				deep++;
//			} else {
//				F[deep + 1] = 0;
//				F[deep]=0;
//				for (int m = 0; m < record[deep + 1].size(); ++m) {
//					int temp = (int) record[deep + 1].get(m);
//					D[temp] = 0;
//				}
//				record[deep + 1].clear();
//				deep--;
//			}
//
//		}
//		motif_degree[i] = motif_num;
//
//	}
//	return motif_degree;
//}