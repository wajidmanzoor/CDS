package hku.algo.cds;

import hku.algo.exist.*;
import hku.algo.prune.*;
import hku.util.*;

import java.util.Date;
import java.util.Map;
import java.util.Queue;

public class ExactTest {

    static String dataset_dir = "./datasets/";
    static String motif_dir = "./motif/";
    static String[] datasets = {"test2"};
    static String[] motifs = { "triangle" };

    public static void main(String[] args) {
        Log stdout = new Log();

        for (String dataset : datasets) {
            try {
                DataReader gReader = new DataReader(dataset_dir + dataset + ".txt", null);
                gReader.readGraph();

                stdout.write("dataset " + dataset + "\n");
                stdout.write("V: "+gReader.graph_size +" E: "+ gReader.count_edge/2+"\n");

                System.out.println(getTime() + " " + dataset);

                for (String motif : motifs) {
                    DataReader mReader = new DataReader(null, motif_dir + motif + ".txt");
                    mReader.readMotif();

                    stdout.write("motif " + motif + "\n");
                    System.out.println(getTime() + " motif " + motif);

                    runAll(stdout, gReader, mReader);
                    mReader = null;
                }


                gReader = null;
                System.gc();

                stdout.write("\n\n\n");

            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void runAll(Log stdout, DataReader g, DataReader m) throws Exception {
        DataReader a = new DataReader(null, null);
        a.Graph = g.Graph;
        a.graph_size = g.gentgraph_size();
        a.Motif = m.Motif;
        a.Motif_Count = m.Motif_Count;
        a.Motif_Type = m.Motif_Type;

        int[][] Graph = a.Graph;
        int[][] Motif = a.Motif;

        long start_time = System.currentTimeMillis();

        CDSdecompose c = new CDSdecompose(Graph, Motif, a.gentgraph_size(),
                                           Motif.length, a.getMotif_Count(), null, null);
        double[][] r_d = c.Decompose();
        long decompose_time = System.currentTimeMillis();
        for (int i = 0; i < r_d.length; i++) {
        for (int j = 0; j < r_d[i].length; j++) {
            System.out.print(r_d[i][j] + " ");
        }
        System.out.println(); // new line after each row
        }





        LocateCore d = new LocateCore(Graph, r_d, a.gentgraph_size());
        DensestCore r_c = d.locate();

        KList b = new KList(r_c.Graph, a.Motif.length);
        b.ListRecord();

        InvalidEdgePruning e = new InvalidEdgePruning(b.Statistic, r_c.Graph, r_c.graph_size);
        e.Prune();

        ComponentDecom f = new ComponentDecom(r_c.Graph, r_c.graph_size, b.Statistic);
        Queue<Component> r_q = f.decompose();

        stdout.write("Num components " + r_q.size() + "\n");

        DensestCore my = new DensestCore(r_c.Graph, r_c.graph_size, 0, 0, 0, 0, r_c.graph_size);
        DynamicExactalgo g_algo = new DynamicExactalgo(r_q, my, Motif[0].length);
        MDS mds = g_algo.DynamicExact();

        stdout.write("Densest Core \n");
        StringBuilder coreStr = new StringBuilder();
        for (int x = 0; x < my.graph_size; x++) {
            coreStr.append(d.reverse_map[x]).append(" ");
        }
        stdout.write(coreStr.toString().trim() + "\n");

        if (mds.s_t_result != null) {
            int[] result = new int[mds.vertex_num];
            int ind = 0;
            for (int x = 0; x < mds.s_t_result.length; x++) {
                int v = mds.s_t_result[x];
                if (v != -1) {
                    if (mds.core.reverse_map != null) {
                        result[ind++] = d.reverse_map[mds.core.reverse_map[v]];
                    } else {
                        result[ind++] = d.reverse_map[v];
                    }
                }
            }

            stdout.write("result \n");
            StringBuilder resultStr = new StringBuilder();
            for (int i = 0; i < ind; i++) {
                resultStr.append(result[i]).append(" ");
            }
            stdout.write(resultStr.toString().trim() + "\n");
        }

        long stop_time = System.currentTimeMillis();

        stdout.write("density " + mds.densest + "\n");
        stdout.write("Num cliques " + mds.motif_num + "\n");
        stdout.write("Num verticies " + mds.vertex_num + "\n");
        stdout.write("Total Time: " + (stop_time - start_time) + "ms\n");
    }

    private static String getTime() {
        return new Date().toLocaleString();
    }
}
