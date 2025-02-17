package hku.algo.maxflow;


import java.util.Arrays;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;

/**
 * compute the Min-Cut, according to the flow networks.
 * 
 * In this class, we used various exiting algorithms. 
 * EK:
 * Preflow: 
 * 
 * @author yukaiqiang
 * @date Nov,23 2017
 */
public class FindMinCut {
	
	/** data structure used to record flow network*/
	public Map<Integer, double[]>[] FlowNetwork=null;
	
	/** array used to save path*/
	public int[] parent=null;
	/** source vertex */
	private int s=0;
	/** sink vertex */
	private int t=0;
	
	/**
	 * 
	 * @param FlowNetwork data structure used to record flow network
	 * @param s source vertex
	 * @param t sink vertex
	 */
	public FindMinCut(Map<Integer,double[]>[] FlowNetwork,int s,int t) {
		this.FlowNetwork=FlowNetwork;
		this.s=s;
		this.t=t;
	}
	
	/**
	 * EK algorithm to solve the max-flow problem
	 * the min-cut(S-T) is saved in the array 'parent[]'
	 * S: parent[i]>-1
	 * T: parent[i]=-1
	 * @return the value of max flow
	 */
	public double EdmondsKarp() {
		parent=new int[FlowNetwork.length];

		//wm: returns the min cut of augmented path and the augmented path 
		double result=augmentPath(parent);

		//wm: store the total flow of network
		double sum=0;

		double temp[];
		
		//wm: break then no augmented path is found
		while(result!=-1) {

			//wm: current vertex start from sink 
			int cur=t;
			double fre=0;

			//wm: backtrack to update the available capacity 
			while(cur!=s) {
				
				//wm: availble capacity before update
				temp=FlowNetwork[parent[cur]].get(cur);

				//wm: 
				fre=temp[0];

				//wm: decreament min cut from the available capacity
				temp[0]=temp[0]-result;

				//wm: update the reverse edge, by increamenting the available capacity 
				// this allows to reverse the flow from the forward edge if needed in future
				temp=FlowNetwork[cur].get(parent[cur]);
				temp[0]=temp[0]+result;

				//wm: set current to the parrent of current
				cur=parent[cur];
				
			}
			

			//wm: increament the total flow of network by min cut of the augmented path.
			sum+=result;
			
			//wm: get new augmented path 
			result=augmentPath(parent);
			
		}
		
		//wm: return the total flow through network
		return sum;
	}
	
	/**
	 * use BFS algorithm to find the augment path. 
	 * @param parent the array to save the path
	 * @return the minimum capacity of each edge in the path
	 */
	private double augmentPath(int[] parent) {

		//wm: saves the minimum cut- or maximum posible flow through the path
		double maxflow=Integer.MAX_VALUE;

		//wm: this stores the parent of the vertex in the flownetwork 
		Arrays.fill(parent, -1);

		//wm: queue for BFS
		Queue<Integer> queue=new LinkedList<Integer>();

		//wm: add source vertex to queue
		queue.add(s);

		//wm: set parent of source as source (as it is the starting point)
		parent[s]=s;
		double temp[];
		
		while(!queue.isEmpty()) {
			int p=queue.poll();

			//wm: if we reach the sink vertex
			if(p==t) {

				//wm: backtracking 
				while(p!=s) {

					//wm: available capacity 
					temp=FlowNetwork[parent[p]].get(p);

					//wm: code to find minimum cut (min available capacity)
					if(maxflow>temp[0])
						maxflow=temp[0];
					//wm: backtrack to parent of vertex in flow network
					p=parent[p];
				}

				//wm: if source is reached break 
				break;
			}
			
			//wm: BFS, iter through childern of vertex in flow network 
			for(Entry<Integer, double[]> entry: FlowNetwork[p].entrySet()) {

				//wm: availible and total capacity
				temp=entry.getValue();

				//wm: if childrn is not visited and has availible capacity 
				if(parent[entry.getKey()]==-1&&temp[0]>0) {

					//wm: mark childern as visited by setting it parent 
					parent[entry.getKey()]=p;

					//wm: add childern to queue
					queue.add(entry.getKey());
					
				}
			}
		}
		
		//wm: if path now found return -1
		if(parent[t]==-1) {
			return -1;
		}

		//wm: if path found return the maximum flow posible (minimum cut) throught the augementated path
		return maxflow;
	}
	
	/**
	 * 
	 * @return the parent
	 */
	public int[] getparent() {

		//wm: give the augmented path 
		return parent;
	}
	
}
