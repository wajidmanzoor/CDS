package hku.algo.exist;

import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Queue;

import hku.algo.prune.Component;
import hku.algo.prune.DensestCore;

public class DynamicExactalgo {
	
	private Queue<Component> queue=null;
	
	private DensestCore Core=null;
	
	private int motif_size;
	
	public DynamicExactalgo(Queue<Component> queue,DensestCore Core,int motif_size) {
		this.queue=queue;
		this.Core=Core;
		this.motif_size=motif_size;
	}
	
	public MDS DynamicExact(){

		//wm: iter through connected components to find the denset one (num motif / component size )
		// index stores the component with higest density
		// low bound stores the higest density
		Component C,index;
		index=queue.iterator().next();
		double low_bound=0;
		Iterator it=queue.iterator();
		while(it.hasNext()) {
			C=(Component) it.next();
			if(low_bound<C.densest) {
				low_bound=C.densest;
				index=C;
			}
		}
		
		//wm: the density of the whole core > the denset component. whoes core density as lower bound
		if(Math.ceil(low_bound)<Math.ceil(Core.densest)) {
			low_bound=Core.densest;
		}
		
		//wm: set upper bound as kmax
		double up_bound=Core.kmax;
		
		Exactalgo exact;
		
		//wm: stores the denset subgraph results
		MDS mds=new MDS(null,index.motif_num,index.graph_size,low_bound);


		//wm: iter over connected components
		while(!queue.isEmpty()) {

			//wm: connected component
			C=queue.remove();
		

			//wm: do not know why this code
			exact=new Exactalgo(C.motif_list, motif_size, C.graph_size, C.motif_degree);
			long test=0;
			for(int mm=0;mm<C.motif_degree.length;++mm) {
				test+=C.motif_degree[mm];
			}

				//wm: create an instance of exact algorithm 
				exact=new Exactalgo(C.motif_list, motif_size, C.Graph.length, C.motif_degree);
				
				//wm: run the exact algorithm of the connected component
				int res[]=exact.Exact(Math.ceil(low_bound), Math.ceil(up_bound)+1, C.motif_num);


				//wm: compute the total motif in the subgraph returned by exact algorithm
				long motif_num=0; int vertex_num=0;
				for(Entry<String,int[]> entry:C.motif_list.entrySet()) {
					int temp[]=entry.getValue();
					int i=0;
					
					//wm: break if any vertex in the motif was removed
					for(;i<temp.length-1;++i) {
						if(res[temp[i]]==-1)
							break;
					}

					//wm: if no vertex in a motif was removed, increament motif num 
					if(i==temp.length-1) {
						motif_num+=temp[i];
					}
						
				}

				//wm: count the verticies in new denset graph
				for(int i=0;i<C.graph_size;++i) {

					//wm: if vertex not removed
					if(res[i]!=-1)
						vertex_num++;
					
				}

				//wm: just to avoid division by zero
				if(vertex_num==0)
					vertex_num=C.graph_size;

				//wm: calculate new motif density density 
				low_bound=(double)motif_num/(vertex_num*1.0);

				//wm: if new motif density is greater than the current one, update the denset graph
				if(low_bound>mds.densest) {
					mds.densest=low_bound;
					mds.motif_num=motif_num;
					mds.s_t_result=res;
					mds.vertex_num=vertex_num;
					mds.core=C;
				}
			}
//		}
		
		return mds;
	}

}
