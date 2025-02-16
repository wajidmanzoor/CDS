package hku.algo.cds;



import java.util.Date;
import java.util.Queue;

import hku.algo.exist.Appalgo;
import hku.algo.exist.DynamicExactalgo;
import hku.algo.exist.Exactalgo;
import hku.algo.exist.MDS;
import hku.algo.findgeneralpattern.FindMotif;
import hku.algo.prune.Component;
import hku.algo.prune.ComponentDecom;
import hku.algo.prune.DensestCore;
import hku.algo.prune.InvalidEdgePruning;
import hku.algo.prune.LocateCore;
import hku.util.Combination;
import hku.util.DataReader;
import hku.util.KCore;
import hku.util.Log;


public class ExactTest {

	//wm: folder path of graph and motif
	static String dataset_doc = "./datasets/";
	static String motif_doc = "./motif/";

	//wm: list of graph and motif filenames 
	static String[] datasets_url = { "graph"};
	static String[] motif_url = { "edge"};

	static int[] motif_d = { 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4 };
	static String result = "./result/AppEfficiency.txt";

	//wm: store adjancey list of graphs
	static DataReader g_data[] = new DataReader[14];

	//wm: store adjancey matrix of motifs
	static DataReader m_data[] = new DataReader[14];
	public static void main(String[] args) {

				
		try {
			//wm: Reads  all graphs in dataset url list
			for (int i = 0; i < datasets_url.length; ++i) {
				
				//wm: stores adjancey list of graphs
				g_data[i] = new DataReader(dataset_doc + datasets_url[i] + ".txt", null);
				g_data[i].readGraph();
			}

			// wm: Reads  all motifs in motif url list
			for (int j = 0; j < motif_url.length; ++j) {

				//wm: sotres adjancey matrix of motif
				m_data[j] = new DataReader(null, motif_doc + motif_url[j] + ".txt");

				m_data[j].readMotif();
			}
			//wm: instance to temporary store the graph data
			DataReader a = new DataReader(null, null);

			//Wm: class that write to extenal file file name set in Config file
			Log stdout = new Log();

						
			 for(int i=0;i<datasets_url.length;++i) {

					stdout.write(datasets_url[i]+"\n");
					Date date = new Date();
					String time = date.toLocaleString();
					System.out.println(time+" "+datasets_url[i]);
					
					date = new Date();
					time = date.toLocaleString();
					
					for(int j=0;j<motif_url.length;++j) {
						stdout.write("  "+motif_url[j]+"\n");
						
						for(int alg=0;alg<2;++alg) {
							if(alg==0) {
								stdout.write("**************CoreExact: Exact+Pruning***************\n");
								System.out.println(time+" **************CoreExact: Exact+Pruning***************"+"dtasets_index: "+i);
								date = new Date();
								time = date.toLocaleString();
								System.out.println(time+" "+motif_url[j]+" motif_index: "+j);
								
								// wm: copy graph adjancey list to another object
								a.Graph=g_data[i].Graph;

								//wm: graph size (num verticies) copied to another object 
								a.graph_size=g_data[i].gentgraph_size();
								
								//wm: copy motif details to another object 
								a.Motif=m_data[j].Motif;
								a.Motif_Count=m_data[j].Motif_Count;
								a.Motif_Type=m_data[j].Motif_Type;
								
								int Graph[][]=a.Graph;
								int Motif[][]=a.Motif;
								
								//wm:: total number of edges
								int counta=0;
								for(int yui=0;yui<Graph.length;++yui){
									counta+=Graph[yui].length;
								}
								
								
								long start_time=System.currentTimeMillis();
								
								//wm: clique core decompose or motif core decompose
								hku.algo.cds.CDSdecompose c=new hku.algo.cds.CDSdecompose(Graph, Motif, a.gentgraph_size(),
										Motif.length,a.getMotif_Count(),null,null);
								
								double r_d[][]=c.Decompose();
//								
								long decompose_time=System.currentTimeMillis();
								
								
								//wm: defines lower bound based of motif degree of maximum average motif density core,
								// deleted vertices with motif degree less than the lower bound
								// created new graph and core for the new subgraph  
								// returns the densest core 
								LocateCore d=new LocateCore(Graph, r_d, a.gentgraph_size());	
								DensestCore r_c=d.locate();
								

								//wm: list all k-cliques and get the updted clique core values
								KList b=new KList(r_c.Graph,a.Motif.length);
								b.ListRecord();
//								

								//wm: remove edges that are not part of any clique
								InvalidEdgePruning e=new InvalidEdgePruning(b.Statistic, r_c.Graph, r_c.graph_size);
								int invalid_edge=e.Prune();
								
								//wm: data structure that savesthe connected component details
								ComponentDecom f=new ComponentDecom(r_c.Graph, r_c.graph_size, b.Statistic);
								
								//wm: data structure that save the denset core based of the lower bound
								DensestCore my=new DensestCore(r_c.Graph,r_c.graph_size,0,0,0,0,r_c.graph_size);
								
								//wm: returns a queue that contains all the connected components (new vertex numbering, motif lists, amd lotif degree)
								Queue<Component> r_q=f.decompose();
								int r_q_size=r_q.size();
							
								
								long pruning_time=System.currentTimeMillis();
								

								
								DynamicExactalgo g=new DynamicExactalgo(r_q,my,Motif[0].length);
								
								MDS mds=g.DynamicExact();
								
								long stop_time=System.currentTimeMillis();
								
								/* info of the CDS/EDS */
								System.out.println("density: "+mds.densest);
								System.out.println("Num of Cliques/edges: "+mds.motif_num);
								System.out.println("Num of Vertices: "+mds.vertex_num);
								stdout.write("      total: "+(stop_time-start_time)+"ms\n");
								//stdout.flush();
							}else if(alg==1) {
								stdout.write("**************Exact (Baseline): Exact without Pruning***************\n");
								System.out.println(time+" **************Exact (Baseline): Exact without Pruning***************"+"dtasets_index: "+i);
								date = new Date();
								time = date.toLocaleString();
								System.out.println(time+" "+motif_url[j]+" motif_index: "+j);
								a.Graph=g_data[i].Graph;
								a.graph_size=g_data[i].gentgraph_size();
								a.Motif=m_data[j].Motif;
								a.Motif_Count=m_data[j].Motif_Count;
								a.Motif_Type=m_data[j].Motif_Type;
								
								int Graph[][]=a.Graph;
								int Motif[][]=a.Motif;
								
								int counta=0;
								for(int yui=0;yui<Graph.length;++yui){
									counta+=Graph[yui].length;
								}
								//System.out.println(Graph.length+"graph_size: "+counta/2);
								
								
								long start_time=System.currentTimeMillis();
								
								KList b=new KList(Graph,Motif.length);
								b.ListRecord();
//								
								int motif_degree[] =new int[Graph.length];
								for(int i1=0;i1<Graph.length;++i1) {
									motif_degree[i1]=(int)b.motif_degree[i1];
								}
								Exactalgo exact=new Exactalgo(b.Statistic, Motif[0].length, Graph.length,motif_degree);
								exact.Exact(0,b.getMotifNum(),b.getMotifNum());
								
								long stop_time=System.currentTimeMillis();
								
								/* info of the CDS/EDS */
								stdout.write("      total: "+(stop_time-start_time)+"ms\n");
								//stdout.flush();
							}
						}
					}
					stdout.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
					
				}
				
//			
				
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	 
	}
	

}
