package hku.algo.cds;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import hku.util.DataReader;
import hku.util.KCore;

public class KList {

    //wm: order of verticis in the sorted vertices array (sorted by core value descending)
    //eg: core [1,2,2,2], sorted vertices =[1,2,3,0], order = [4,1,2,3] 
    private int[] order;
    
    // wm: store neighbors of each vertex that are after the vertex in order array
    //order = [[],[0,2,3],[3],[]]
    //wm: stored DAG representation as an adjacency list
    private int[][] graph;

    //wm: clique number or motif size
    private int k;
    
    //wm: graph size or num of verticies
    private static int graph_size;

    //wm: 
    private int degree[];
    private int label[];
    public int motif_num=0;
    
    //wm : Clique degree
    public long[] motif_degree;

    //wm: Adjancy list of graph
    private int[][] GenGraph;

    //wm: stores all cliques string is clique representation and then array is all the cliques
    public Map<String, int[]> Statistic=new HashMap<String, int[]>();
    
    public ArrayList<String> Adjlist[];
    
    
    public KList(int[][] graph,int k) {
            this.GenGraph=graph;
            this.k=k;
            
            graph_size=graph.length;
            this.graph=new int[graph_size][];
            order=new int[graph_size];
            degree=new int[graph_size];
            label=new int[graph_size];
            Arrays.fill(label, k);


            motif_degree=new long[graph_size];
            Arrays.fill(motif_degree, 0);
    }
    
    public void getListingOrder() {
            KCore b=new KCore(GenGraph);
            //wm: simple core decompose algorithm
            b.decompose();
            
            // wm: verticies sorted by core value, descending order
            int temp_arr[]=b.obtainReverseCoreArr();

            // wm : give the order of each vertex by core value (descending ), from 1 to n
            for(int i=0;i<graph_size;++i) {
                    order[temp_arr[i]]=i+1;
            }
            
    }
    
    public void GenerateDAG() {
        /**
         *  Converts Undirected Graph into Directed Acyclic Graph
         *  Edges only go forward (i.e., from a lower-ordered vertex to a higher-ordered one).
         *  Reducing search space by ensuring we only explore forward-moving edges in the order.
        **/
            
            int count;
            for(int i=0;i<graph_size;++i) {
                    count=0;
                    //wm: count number of neighbors of vertex that are after the vertex in order array
                    for(int j=0;j<GenGraph[i].length;++j) {
                            if(order[i]<order[GenGraph[i][j]]) {
                                    count++;
                            }
                    }
                    int arr[]=new int[count];
                    // wm: set that count to degree
                    degree[i]=count;
                    count=0;
                    //wm: get the neighbors of vertex that are after the vertex in order array
                    for(int j=0;j<GenGraph[i].length;++j) {
                            if(order[i]<order[GenGraph[i][j]]) {
                                    arr[count]=GenGraph[i][j];
                                    count++;
                            }
                    }
                    graph[i]=arr;
            }
            
    }
    
    public void ListingRecord(int k,ArrayList<Integer> c,ArrayList<Integer> arr) {

            //wm: lists all k cliques
            if(k==2) {
                    String a="";
                    //wm : get the string for hash map 
                    for(int m=0;m<c.size();++m) {
                            a+=c.get(m)+" ";
                    }

                 //wm: Tracks the number of motif instances (edges).
                    int multi=0;

                    //wm: iter over all candidate verticies
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);

                            //wm: iter over valid neighbors of candidate vertex 
                            for(int j=0;j<degree[temp];++j) {
                                    
                                    //wm: add to hasp map string 
                                    a=a+temp+" "+graph[temp][j];

                                    //wm: increament number of instances
                                    multi++;

                                    //wm: increament gloabl number of motif
                                    motif_num++;

                                    //wm: increase motif degree of neighbor
                                    motif_degree[graph[temp][j]]++;

                                    //wm: increase motif degree of candidate vertex
                                    motif_degree[temp]++;

                                    //wm: store the k-clique verticies
                                    int temp_arr[]=new int[this.k+1];

                                    //WM: get all the partial clique verticies 
                                    for(int m=0;m<c.size();++m) {
                                            temp_arr[m]=c.get(m);
                                    }

                                    //wm: add candidate vertex and its valid neighbor
                                    temp_arr[this.k-2]=temp;
                                    temp_arr[this.k-1]=graph[temp][j];

                                    //wm: mark as completed clique
                                    temp_arr[this.k]=1;
                                    
                                    //wm: add to map
                                    this.Statistic.put(a, temp_arr);
                            }
                    }
                
                    //wm : increase motif degree of all verticies in C
                    for(int m=0;m<c.size();++m) {
                            int temp=c.get(m);
                            motif_degree[temp]+=multi;
                    }
                    
            }else {
                    
                    //wm: iter through all verticies in candidate set 
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);
                            
                            //wm: stores DAG neighbors of candidate set that are valid. 
                            ArrayList<Integer> arr_n=new ArrayList<Integer>();
                            for(int j=0;j<graph[temp].length;++j) {
                                    if(label[graph[temp][j]]==k) {

                                            //wm: maek them as processed
                                            label[graph[temp][j]]=k-1;
                                            
                                            //wm: add to ar_n
                                            arr_n.add(graph[temp][j]);
                                    }                                               
                            }
                            
                            //wm: iter through valid neghbors 
                            for(int j=0;j<arr_n.size();++j) {
                                    int arr_temp=arr_n.get(j);

                                    int index=0;

                                    //wm: reorder the neighors of valid neighbors to put valid neighbor first

                                    //wm: iter through valid neighborsfrom end
                                    for(int m=graph[arr_temp].length-1;m>index;--m) {
                                            if(label[graph[arr_temp][m]]==k-1) {

                                                    //wm: get index of first non valid neighbors ( which is before the index of valid neighbor)
                                                    while(index<m&&label[graph[arr_temp][index]]==k-1) {
                                                            index++;
                                                    }
                                                    if(label[graph[arr_temp][index]]!=k-1) {
                                                            //wm: swap the valid and non valid neighbor
                                                            int temp1=graph[arr_temp][m];
                                                            graph[arr_temp][m]=graph[arr_temp][index];
                                                            graph[arr_temp][index]=temp1;
                                                            
                                                    }
                                            }
                                    }

                                    //wm: coz indexing starts from zero and count will be +1 of index
                                    if(graph[arr_temp].length!=0&&label[graph[arr_temp][index]]==k-1)
                                            index++;
                                    //wm: store count of valid neighbors
                                    degree[arr_temp]=index;
                            }
                            
                            //wm: add vertex fron candidate set to  partial clique set 
                            c.add(arr.get(i));

                            //wm: recursively call function for k-1 clique
                            ListingRecord(k-1,c,arr_n);

                            //wm: remove vertex from c and update the label for next vertex
                            c.remove(arr.get(i));
                            for(int j=0;j<arr_n.size();++j) {
                                    int arr_temp=arr_n.get(j);
                                    label[arr_temp]=k;
                            }
                    }
            }
            
    }
    
    
    
    public void ListingAdj(int k,ArrayList<Integer> c,ArrayList<Integer> arr) {

        //wm: lists all k cliques
        if(k==2) {
                //System.out.println(">>>");

                //wm : get the string for hash map 
                String a="";
                for(int m=0;m<c.size();++m) {
                        a+=c.get(m)+" ";
                }

                 //wm: Tracks the number of motif instances (edges).
                int multi=0;

               //wm: iter over all candidate verticies
                for(int i=0;i<arr.size();++i) {
                        int temp=arr.get(i);

                        //wm: iter over valid neighbors of candidate vertex 
                        for(int j=0;j<degree[temp];++j) {
                               
                               //wm: add to hasp map string 
                                String cc=a+temp+" "+graph[temp][j];

                               //wm: increament number of instances
                                multi++;

                               //wm: increament gloabl number of motif
                                motif_num++;

                               //wm: increase motif degree of neighbor
                                motif_degree[graph[temp][j]]++;

                               //wm: increase motif degree of candidate vertex
                                motif_degree[temp]++;

                                //wm: store the k-clique verticies
                                int temp_arr[]=new int[this.k+1];

                                //WM: get all the partial clique verticies 
                                for(int m=0;m<c.size();++m) {
                                        temp_arr[m]=c.get(m);
                                }

                                //wm: add candidate vertex and its valid neighbor
                                temp_arr[this.k-2]=temp;
                                temp_arr[this.k-1]=graph[temp][j];

                                //wm: mark as completed clique
                                temp_arr[this.k]=1;
                                //wm: add to map
                                this.Statistic.put(cc, temp_arr);

                                //wm: store all cliques keys in a adjancey list 
                                for(int m=0;m<this.k;++m) {
                                		Adjlist[temp_arr[m]].add(cc);
                                }
                        }
                }

                //wm : increase motif degree of all verticies in C
                for(int m=0;m<c.size();++m) {
                        int temp=c.get(m);
                        motif_degree[temp]+=multi;
                }
                
        }else {
                
                //wm: iter through all verticies in candidate set 
                for(int i=0;i<arr.size();++i) {
                        int temp=arr.get(i);
                        
                        //wm: stores DAG neighbors of candidate set that are valid. 
                        ArrayList<Integer> arr_n=new ArrayList<Integer>();
                        for(int j=0;j<graph[temp].length;++j) {
                                
                                if(label[graph[temp][j]]==k) {
                                        
                                        //wm: mark them as processed
                                        label[graph[temp][j]]=k-1;
                                        
                                        //wm: add to ar_n
                                        arr_n.add(graph[temp][j]);
                                }                                               
                        }
                        
                        //wm: iter through valid neghbors 
                        for(int j=0;j<arr_n.size();++j) {

                                int arr_temp=arr_n.get(j);
                                int index=0;

                                //wm: reorder the neighors to put valid neighbor first

                                //wm: iter through valid neighborsfrom end
                                for(int m=graph[arr_temp].length-1;m>index;--m) {
                                        if(label[graph[arr_temp][m]]==k-1) {

                                                //wm: get index of first non valid neighbors ( which is before the index of valid neighbor)
                                                while(index<m&&label[graph[arr_temp][index]]==k-1) {
                                                        index++;
                                                }
                                                if(label[graph[arr_temp][index]]!=k-1) {
                                                        
                                                        //wm: swap the valid and non valid neighbor
                                                        int temp1=graph[arr_temp][m];
                                                        graph[arr_temp][m]=graph[arr_temp][index];
                                                        graph[arr_temp][index]=temp1;
                                                        
                                                }
                                        }
                                }

                                //wm: coz indexing starts from zero and count will be +1 of index
                                if(graph[arr_temp].length!=0&&label[graph[arr_temp][index]]==k-1)
                                        index++;
                                //wm: store count of valid neighbors
                                degree[arr_temp]=index;
                        }
                        
                        //wm: add vertex fron candidate set to  partial clique set
                        c.add(arr.get(i));

                        //wm: recursively call function for k-1 clique
                        ListingAdj(k-1,c,arr_n);

                        //wm: remove vertex from c and update the label for next candidate vertex
                        c.remove(arr.get(i));
                        for(int j=0;j<arr_n.size();++j) {
                                int arr_temp=arr_n.get(j);
                                label[arr_temp]=k;
                        }
                }
        }
        
}
    
    
    public void Listing(int k,ArrayList<Integer> c,ArrayList<Integer> arr) {

        // Returns all k- cliques

            //wm: 2 -clique (edge)
            if(k==2) {
                    
                    String a="";
                    //wm : get the string for hash map 
                    for(int m=0;m<c.size();++m) {
                            a+=c.get(m)+" ";
                    }

                    //wm: Tracks the number of motif instances (edges).
                    int multi=0;

                    //wm: iter over all candidate verticies
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);
                            //wm: iter over valid neighbors of candidate vertex 
                            for(int j=0;j<degree[temp];++j) {

                                    //wm: increament number of instances
                                    multi++;
                                    //wm: increament gloabl number of motif
                                    motif_num++;

                                    //wm: increase motif degree of neighbor
                                    motif_degree[graph[temp][j]]++;

                                    //wm: increase motif degree of candidate vertex
                                    motif_degree[temp]++;
                            }
                    }

                    //wm : increase motif degree of all verticies in C
                    for(int m=0;m<c.size();++m) {
                            int temp=c.get(m);
                            motif_degree[temp]+=multi;
                    }
            //wm: rest
            }else {
                    

                    //wm: iter through all verticies in candidate set 
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);
                            
                            //wm: stores DAG neighbors of candidate set that are valid.
                            ArrayList<Integer> arr_n=new ArrayList<Integer>();
                            for(int j=0;j<graph[temp].length;++j) {
                                    if(label[graph[temp][j]]==k) {

                                            //wm: maek them as processed
                                            label[graph[temp][j]]=k-1;
                                            
                                            //wm: add to ar_n
                                            arr_n.add(graph[temp][j]);
                                    }                                               
                            }
                            
                            //wm: iter through valid neghbors 
                            for(int j=0;j<arr_n.size();++j) {
                                    int arr_temp=arr_n.get(j);
                                    int index=0;

                                    //wm: reorder the neighors to put valid neighbor first

                                    //wm: iter through valid neighborsfrom end
                                    for(int m=graph[arr_temp].length-1;m>index;--m) {
                                            if(label[graph[arr_temp][m]]==k-1) {

                                                    //wm: get index of first non valid neighbors ( which is before the index of valid neighbor)
                                                    while(index<m&&label[graph[arr_temp][index]]==k-1) {
                                                            index++;
                                                    }
                                                    if(label[graph[arr_temp][index]]!=k-1) {

                                                            //wm: swap the valid and non valid neighbor
                                                            int temp1=graph[arr_temp][m];
                                                            graph[arr_temp][m]=graph[arr_temp][index];
                                                            graph[arr_temp][index]=temp1;
                                                            
                                                    }
                                            }
                                    }

                                    //wm: coz indexing starts from zero and count will be +1 of index
                                    if(graph[arr_temp].length!=0&&label[graph[arr_temp][index]]==k-1)
                                            index++;

                                    //wm: store count of valid neighbors
                                    degree[arr_temp]=index;
                            }
                            
                            //wm: add vertex fron candidate set to  partial clique set 
                            c.add(arr.get(i));

                            //wm: recursively call function for k-1 clique
                            Listing(k-1,c,arr_n);

                            //wm: remove vertex from c and update the label for next vertex
                            c.remove(arr.get(i));
                            for(int j=0;j<arr_n.size();++j) {
                                    int arr_temp=arr_n.get(j);
                                    label[arr_temp]=k;
                            }
                    }
            }
            
    }
    
    
    
    public void Listing(int k,ArrayList<Integer> c,ArrayList<Integer> arr,int map) {

            //list all k-cliques that include a desrired vertex
            if(k==2) {
                    boolean onenode=false;
                    String a="";

                    //wm: set one node to true if partial clique include the desired vertex
                    for(int m=0;m<c.size();++m) {
                            a+=c.get(m)+" ";
                            if(c.get(m)==map) {
                                    onenode=true;
                            }
                    }

                    //wm: Tracks the number of motif instances (edges).
                    int multi=0;

                    //wm: iter over all candidate verticies
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);

                            //wm: iter over valid neighbors of candidate vertex 
                            for(int j=0;j<degree[temp];++j) {
//                      
                                    //wm: if cliques containes the desired vertex
                                    if(onenode||temp==map||graph[temp][j]==map) {

                                            //wm: increament number of instances
                                            multi++;

                                            //wm: increament gloabl number of motif
                                            motif_num++;

                                            //wm: increase motif degree of neighbor
                                            motif_degree[graph[temp][j]]++;

                                            //wm: increase motif degree of candidate vertex
                                            motif_degree[temp]++;
                                    }
                                    
                            }
                    }

                    //wm : increase motif degree of all verticies in C
                    for(int m=0;m<c.size();++m) {
                            int temp=c.get(m);
                            motif_degree[temp]+=multi;
                    }
                    
            }else {
                    
                    //wm: iter through all verticies in candidate set 
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);
                            
                            //wm: stores DAG neighbors of candidate set that are valid.
                            ArrayList<Integer> arr_n=new ArrayList<Integer>();
                            for(int j=0;j<graph[temp].length;++j) {
                                    if(label[graph[temp][j]]==k) {

                                            //wm: maek them as processed
                                            label[graph[temp][j]]=k-1;
                                            
                                            //wm: add to ar_n
                                            arr_n.add(graph[temp][j]);
                                    }                                               
                            }
                            
                            //wm: iter through valid neghbors 
                            for(int j=0;j<arr_n.size();++j) {

                                    int arr_temp=arr_n.get(j);
                                    int index=0;

                                    //wm: reorder the neighors to put valid neighbor first

                                    //wm: iter through valid neighborsfrom end
                                    for(int m=graph[arr_temp].length-1;m>index;--m) {
                                            if(label[graph[arr_temp][m]]==k-1) {

                                                    //wm: get index of first non valid neighbors ( which is before the index of valid neighbor)
                                                    while(index<m&&label[graph[arr_temp][index]]==k-1) {
                                                            index++;
                                                    }
                                                    if(label[graph[arr_temp][index]]!=k-1) {

                                                            //wm: swap the valid and non valid neighbor
                                                            int temp1=graph[arr_temp][m];
                                                            graph[arr_temp][m]=graph[arr_temp][index];
                                                            graph[arr_temp][index]=temp1;
                                                            
                                                    }
                                            }
                                    }

                                    //wm: coz indexing starts from zero and count will be +1 of index
                                    if(graph[arr_temp].length!=0&&label[graph[arr_temp][index]]==k-1)
                                            index++;
                                    //wm: store count of valid neighbors
                                    degree[arr_temp]=index;
                            }
                            
                            //wm: add vertex fron candidate set to  partial clique set
                            c.add(arr.get(i));

                            //wm: recursively call function for k-1 clique
                            Listing(k-1,c,arr_n,map);

                            //wm: remove vertex from c and update the label for next vertex
                            c.remove(arr.get(i));
                            for(int j=0;j<arr_n.size();++j) {
                                    int arr_temp=arr_n.get(j);
                                    label[arr_temp]=k;
                            }
                    }
            }
            
    }
    
    
    public void Listing(int k,ArrayList<Integer> c,ArrayList<Integer> arr,int map[]) {

            if(k==2) {
                    boolean onenode=false;
                    String a="";
                    for(int m=0;m<c.size();++m) {
                            a+=c.get(m)+" ";
                            if(map[c.get(m)]==1) {
                                    onenode=true;
                            }
                    }
                    int multi=0;
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);
                            for(int j=0;j<degree[temp];++j) {
//                                  System.out.println(a+temp+" "+graph[temp][j]);
                                    if(onenode||map[temp]==1||map[graph[temp][j]]==1) {
                                            multi++;
                                            motif_num++;
                                            motif_degree[graph[temp][j]]++;
                                            motif_degree[temp]++;
                                    }
                                    
                            }
                    }
                    for(int m=0;m<c.size();++m) {
                            int temp=c.get(m);
                            motif_degree[temp]+=multi;
                    }
                    
            }else {
                    
                    for(int i=0;i<arr.size();++i) {
                            int temp=arr.get(i);
                            //int count=0;
                            
                            ArrayList<Integer> arr_n=new ArrayList<Integer>();
                            for(int j=0;j<graph[temp].length;++j) {
                                    //System.out.println("****"+graph[temp][j]+" "+" "+label[graph[temp][j]]+" "+k);
                                    if(label[graph[temp][j]]==k) {
                                            label[graph[temp][j]]=k-1;
                                            
                                            //count++;
                                            arr_n.add(graph[temp][j]);
                                    }                                               
                            }
                            
                            for(int j=0;j<arr_n.size();++j) {
                                    //int count=0;
                                    int arr_temp=arr_n.get(j);
                                    int index=0;
                                    for(int m=graph[arr_temp].length-1;m>index;--m) {
                                            if(label[graph[arr_temp][m]]==k-1) {
                                                    while(index<m&&label[graph[arr_temp][index]]==k-1) {
                                                            index++;
                                                    }
                                                    if(label[graph[arr_temp][index]]!=k-1) {
                                                            int temp1=graph[arr_temp][m];
                                                            graph[arr_temp][m]=graph[arr_temp][index];
                                                            graph[arr_temp][index]=temp1;
                                                            
                                                    }
                                            }
                                    }
                                    if(graph[arr_temp].length!=0&&label[graph[arr_temp][index]]==k-1)
                                            index++;
                                    degree[arr_temp]=index;
                            }
                            
                            c.add(arr.get(i));
                            Listing(k-1,c,arr_n,map);
                            c.remove(arr.get(i));
                            for(int j=0;j<arr_n.size();++j) {
                                    int arr_temp=arr_n.get(j);
                                    label[arr_temp]=k;
                            }
                    }
            }
            
    }
    
    public void ListFast() {
        
            //wm : apply core decompose to get vertices sorted by core values (descending)
            // also order, order of each vertex in sorted array of vertices by core values (descending)
            getListingOrder();

            //wm: get number of neighbors of each vertex that are after the vertex in order array (i.e who core value is less that vertex core value)
            //wm: the count of those neighbors will be the new degree
            GenerateDAG();

            //wm: store verticies in the partial clique being formed
            ArrayList<Integer> c=new ArrayList<Integer>();

            //wm : The candidate vertices that can be added to extend the clique.
            // Initially includes all verticies
            ArrayList<Integer> arr=new ArrayList<Integer>();                
            for(int i=0;i<graph_size;++i) {
                    arr.add(i);
            }

            //wm: list all the K-cliques 
            Listing(k, c, arr);
    }
    
    
    public void ListRecord() {

            //wm : apply core decompose to get vertices sorted by core values (descending)
            // also order, order of each vertex in sorted array of vertices by core values (descending)
            getListingOrder();


            //wm: get number of neighbors of each vertex that are after the vertex in order array (i.e who core value is less that vertex core value)
            //wm: the count of those neighbors will be the new degree
            GenerateDAG();
            ArrayList<Integer> c=new ArrayList<Integer>();
            ArrayList<Integer> arr=new ArrayList<Integer>();                
            for(int i=0;i<graph_size;++i) {
                    arr.add(i);
            }
            ListingRecord(k, c, arr);
    }
    
    public void ListAdj() {
        getListingOrder();
        GenerateDAG();
        ArrayList<Integer> c=new ArrayList<Integer>();
        ArrayList<Integer> arr=new ArrayList<Integer>();                
        for(int i=0;i<graph_size;++i) {
                arr.add(i);
        }
        Adjlist=new ArrayList[graph_size];
        for(int i=0;i<graph_size;++i) {
        		Adjlist[i]=new ArrayList<String>();
        }
        ListingAdj(k, c, arr);
}
    
    
    public void ListOne(int a) {

            //wm : apply core decompose to get vertices sorted by core values (descending)
            // also order, order of each vertex in sorted array of vertices by core values (descending)
            getListingOrder();

            //wm: get number of neighbors of each vertex that are after the vertex in order array (i.e who core value is less that vertex core value)
            //wm: the count of those neighbors will be the new degree
            GenerateDAG();
            
            // wm: partially cliques 
            ArrayList<Integer> c=new ArrayList<Integer>();

            //wm: candidate verticies
            ArrayList<Integer> arr=new ArrayList<Integer>();                
            for(int i=0;i<graph_size;++i) {
                    arr.add(i);
            }
            Listing(k, c, arr,a);
    }
    
    public void Listbash(int a[]) {
            getListingOrder();
            GenerateDAG();
            
            ArrayList<Integer> c=new ArrayList<Integer>();
            ArrayList<Integer> arr=new ArrayList<Integer>();                
            for(int i=0;i<graph_size;++i) {
                    arr.add(i);
            }
            Listing(k, c, arr,a);
    }
    
    
    
    
    public int getMotifNum() {
            return this.motif_num;
    }
    public long[] getMotifDegree() {
            return this.motif_degree;
    }
    
	
	public static void main(String[] args) {
		
		DataReader a=new DataReader("./datasets/ca-con.txt","./motif/edge.txt");
		//DataReader a=new DataReader("./motif/3triangle.txt","./motif/3triangle.txt");
		int Graph[][]=a.readGraph();
		int Motif[][]=a.readMotif();
		
		KList k=new KList(Graph, 3);
//		k.getListingOrder();
//		k.GenerateDAG();
//		for(int i=0;i<graph_size;++i) {
//			System.out.print(k.order[i]+" ");
//		}
//		System.out.println();
//		System.out.println("*************");
//		for(int i=0;i<k.graph_size;++i) {
//			System.out.print(i+" ");
//			for(int j=0;j<k.graph[i].length;++j) {
//				System.out.print(k.graph[i][j]+" ");
//			}
//			System.out.println();
//		}
//		System.out.println("******");
//		System.out.println("*************");
//		ArrayList<Integer> c=new ArrayList<Integer>();
//		ArrayList<Integer> arr=new ArrayList<Integer>();
//		
//		for(int i=0;i<graph_size;++i) {
//			arr.add(i);
//		}
//		k.Listing(4, c, arr);
		//k.ListFast();
		//k.ListOne(33);
		//k.ListOne(1);
		int []aa= new int[a.graph_size];
		Arrays.fill(aa, 0);
		aa[33]=1;
		aa[855]=1;
		k.Listbash(aa);
		//System.out.println(k.getMotifDegree()[33]);
		
	}

}
