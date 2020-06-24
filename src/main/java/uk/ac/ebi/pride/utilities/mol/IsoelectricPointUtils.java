package uk.ac.ebi.pride.utilities.mol;

import java.util.HashMap;
import java.util.Map;

/**
 * IsoelectricPointUtils is used to calculate the theoretical isoelectric point of a peptide
 * <p/>
 * At the moment we only support one method, and the peptide PTMs are not taken into account
 *
 * <p/>
 * <p/>
 * @author ypriverol
 * @author rwang
 * Date: 03/08/2011
 * Time: 13:55
 */
public class IsoelectricPointUtils {

    private final static BjellpI bjellpI = new BjellpI();

    public static double calculate(String peptideSeq) {
        peptideSeq = peptideSeq.replace("*","");
        peptideSeq = replaceSpecialAA(peptideSeq);
        return bjellpI.calculate(peptideSeq);
    }

    public static String replaceSpecialAA(String seq){
        for(int i=0; i < seq.length(); i++){
            if(!bjellpI.Cterm_pI_expasy.containsKey(String.valueOf(seq.charAt(i)))){
                seq = seq.replaceFirst(String.valueOf(seq.charAt(i)),"");
                i--;
            }
        }

        return seq;
    }

    public static class BjellpI {
        private Map<String,Double> Cterm_pI_expasy = new HashMap<>();
        private Map<String,Double> Nterm_pI_expasy = new HashMap<>();
        private Map<String, Double> sideGroup_pI_expasy = new HashMap<>();
        private double FoRmU = 0.0D;

        public BjellpI(){
            fillMaps();
        }

        public double calculate(String sequence){
            final double epsilon = 0.001;
	        final int iterationMax = 10000;
            int counter = 0;
            double pHs = -2;
            double pHe = 16;
            double pHm;

            while ((counter < iterationMax) && (Math.abs(pHs - pHe) >= epsilon)) {
                pHm = (pHs + pHe) / 2;
                //System.out.println("[" + pHs + ", " + pHm + "]");
                final double pcs = getpI(sequence, Nterm_pI_expasy, Cterm_pI_expasy, sideGroup_pI_expasy, pHs);
		        final double pcm = getpI(sequence, Nterm_pI_expasy, Cterm_pI_expasy, sideGroup_pI_expasy, pHm);
                if (pcs < 0) {
                    return pHs;
                }
                if (((pcs < 0) && (pcm > 0)) || ((pcs > 0) && (pcm < 0))) {
                    pHe = pHm;
                } else {
                    pHs = pHm;
                }
                counter++;
            }
            //System.out.println("[" + pHs + "," + pHe + "], iteration = " + counter);
            double pHround = Math.round(((pHs + pHe) / 2) * 100.0D);
	        return (pHround / 100.0D);
        }

        private double getpI(String sequence, Map<String,Double> AApI_n, Map<String,Double> AApI_c, Map<String,Double> AApI_side, double PH){
            String sideAA;
            double pHpK;
            double FoRmU = 0.0D;

            String ntermAA = String.valueOf(sequence.charAt(0));
            pHpK = PH - Double.valueOf(AApI_n.get(ntermAA).toString());
            FoRmU += 1.0D / (1.0D + Math.pow(10.0D, pHpK));
            String cterm = String.valueOf(sequence.charAt(sequence.length() - 1));
            pHpK = Double.valueOf(AApI_c.get(cterm).toString()) - PH;
            FoRmU += -1.0D / (1.0D + Math.pow(10.0D, pHpK));

            for (int t = 0; t < sequence.length(); ++t) {
                sideAA = String.valueOf(sequence.charAt(t));
                if (AApI_side.containsKey(sideAA)) {
                    double value = Double.valueOf(AApI_side.get(sideAA).toString());
                    if (value < 0.0D) {
                        pHpK = PH + value;
                        FoRmU += 1.0D / (1.0D + Math.pow(10.0D, pHpK));
                    } else {
                        pHpK = value - PH;
                        FoRmU += -1.0D / (1.0D + Math.pow(10.0D, pHpK));
                    }
                }
            }
            return FoRmU;
        }

        private void fillMaps() {
            this.Cterm_pI_expasy.put("A", 3.55D);
            this.Cterm_pI_expasy.put("R", 3.55D);
            this.Cterm_pI_expasy.put("N", 3.55D);
            this.Cterm_pI_expasy.put("D", 4.55D);
            this.Cterm_pI_expasy.put("C", 3.55D);
            this.Cterm_pI_expasy.put("E", 4.75D);
            this.Cterm_pI_expasy.put("Q", 3.55D);
            this.Cterm_pI_expasy.put("G", 3.55D);
            this.Cterm_pI_expasy.put("H", 3.55D);
            this.Cterm_pI_expasy.put("I", 3.55D);
            this.Cterm_pI_expasy.put("L", 3.55D);
            this.Cterm_pI_expasy.put("K", 3.55D);
            this.Cterm_pI_expasy.put("M", 3.55D);
            this.Cterm_pI_expasy.put("F", 3.55D);
            this.Cterm_pI_expasy.put("P", 3.55D);
            this.Cterm_pI_expasy.put("S", 3.55D);
            this.Cterm_pI_expasy.put("T", 3.55D);
            this.Cterm_pI_expasy.put("W", 3.55D);
            this.Cterm_pI_expasy.put("Y", 3.55D);
            this.Cterm_pI_expasy.put("V", 3.55D);
            this.Nterm_pI_expasy.put("A", 7.59D);

            this.Nterm_pI_expasy.put("R", 7.5D);
            this.Nterm_pI_expasy.put("N", 7.5D);
            this.Nterm_pI_expasy.put("D", 7.5D);
            this.Nterm_pI_expasy.put("C", 7.5D);
            this.Nterm_pI_expasy.put("E", 7.7D);
            this.Nterm_pI_expasy.put("Q", 7.5D);
            this.Nterm_pI_expasy.put("G", 7.5D);
            this.Nterm_pI_expasy.put("H", 7.5D);
            this.Nterm_pI_expasy.put("I", 7.5D);
            this.Nterm_pI_expasy.put("L", 7.5D);
            this.Nterm_pI_expasy.put("K", 7.5D);
            this.Nterm_pI_expasy.put("M", 7.0D);
            this.Nterm_pI_expasy.put("F", 7.5D);
            this.Nterm_pI_expasy.put("P", 8.359999999999999D);
            this.Nterm_pI_expasy.put("S", 6.93D);
            this.Nterm_pI_expasy.put("T", 6.82D);
            this.Nterm_pI_expasy.put("W", 7.5D);
            this.Nterm_pI_expasy.put("Y", 7.5D);
            this.Nterm_pI_expasy.put("V", 7.44D);

            this.sideGroup_pI_expasy.put("R", -12.0D);
            this.sideGroup_pI_expasy.put("D", 4.05D);
            this.sideGroup_pI_expasy.put("C", 9.0D);
            this.sideGroup_pI_expasy.put("E", 4.45D);
            this.sideGroup_pI_expasy.put("H", -5.98D);
            this.sideGroup_pI_expasy.put("K", -10.0D);
            this.sideGroup_pI_expasy.put("Y", 10.0D);
  }
    }
}
