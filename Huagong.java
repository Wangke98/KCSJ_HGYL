import java.util.HashMap;
import java.util.Map;

public class HuaGong {

    private static final int zww = 72;//正戊烷摩尔质量
    private static final int zjw = 86;//正己烷摩尔质量
    private static final double product_concentration = 0.95;//产品浓度
    private static final int work_day = 330;//工作日
    private static final int work_hours = 24;//工作时
    private static final double reflux_ratio = 1.4;//回流比
    private static final double P0 = 101.325;//大气压强
    private static final double Pg = 5;//表压
    private static final double K = 237.15;//绝对零度

    int year_process_mount = 45000;//年处理量
    double mass_frac = 0.6;//进料组成

    double F;//原料液流量
    double D;//塔顶产品流量
    double W;//塔底残液流量
    double Xf;//原料液的摩尔分率
    double Xd;//塔顶产品的摩尔分率
    double Xw;//塔底产品的摩尔分率
    double a;//相对挥发度
    double[] Xn;//存放X1,X2,X3...Xn
    double[] Yn;//存放Y1,Y2,Y3...Yn

    double[] temp;//温度数组
    double[] array1;//正戊烷
    double[] array2;//正己烷


    void chapter_2_1() {

        System.out.println("**********2.1的计算开始**********");

        calculate_mole_frac();
        double Mf = calculate_average_molar_mass(Xf);
        double Md = calculate_average_molar_mass(Xd);
        double Mw = calculate_average_molar_mass(Xw);
        material_balance(Mf, Xw, Xf, Xd);

        System.out.println("Xf=" + Xf + "\tXd=" + Xd + "\tXw=" + Xw);
        System.out.println("Mf=" + Mf + "\tMd=" + Md + "\tMw=" + Mw);
        System.out.println("W=" + W + "\tD=" + D + "\tF=" + F);
    }

    Map<String, Double> chapter_2_2() {

        System.out.println("**********2.2的计算开始**********");

        double tF, tD, tW, t1, t2, x1, y1, x2, y2;

        tF = 45 + (Xf - 0.62) * (50 - 45) / (0.45 - 0.62);
        tD = 40 + (Xd - 0.82) * (40 - 36.1) / (0.82 - 1);
        tW = 68.7 + (Xw - 0) * (68.7 - 65) / (0 - 0.07);

        t1 = (tF + tD) / 2;
        t2 = (tF + tW) / 2;

        x1 = quick_calculate(5, 0.62, 0.82, t1 - 40);
        y1 = quick_calculate(5,0.83,0.93,t1 - 40);
        x2 = quick_calculate(5,0.18,0.31,t2 - 55);
        y2 = quick_calculate(5,0.38,0.57,t2 - 55);

        HashMap<String, Double> map = new HashMap<>();
        map.put("tF", tF);
        map.put("tD", tD);
        map.put("tW", tW);
        map.put("t1", t1);
        map.put("t2", t2);
        map.put("x1", x1);
        map.put("y1", y1);
        map.put("x2", x2);
        map.put("y2", y2);
        return map;

    }

    Map<String, Double> chapter_3_1(double t1, double t2) {
        System.out.println("**********3.1的计算开始**********");
        double yP, Xp, Rm, R, L, V, L1, V1, a1, b1, a2, b2;//a1,a2斜率，b1,b2截距
        int tbs = 0, djk = 0, llb = 0;//塔板数，第几块为进料板，共有几块理论板

        calculate_relative_volatility(t1, t2);

        Xp = Xf;
        yP = (a * Xp) / (1 + (a - 1) * Xp);
        Rm = (Xd - yP) / (yP - Xp);
        R = reflux_ratio * Rm;
        L = R * D;
        V = (R + 1) * D;
        L1 = L + F;
        V1 = V;

        a1 = L / V;
        b1 = D / V * Xd;

        a2 = L1 / V1;
        b2 = W / V1 * Xw;

        System.out.println("y=" + a1 + "x+" + b1);
        System.out.println("y'=" + a2 + "x-" + b2);
        System.out.println("x=y/(" + a + "-(" + a + "-1)*y)");

        //y = a1 * x + b1 | y' = a2 * x - b2 | x = y / (a - (a -1) * y)

        double[][] num = new double[20][20]; //设计一个二维数组num[]来装Xn与Yn | Xn相当于num[n][0] | Yn相当于num[0][n]
        //先将数组中的数据赋值很大，以便满足后面的条件
        for (int i = 0; i < 20; i++) {
            for (int j = 0; j < 20; j++) {
                if (i == 0 || j == 0) {
                    num[i][j] = 1;
                }
            }
        }
        double x;
        double y = Xd;
        num[0][1] = y;
        int count = 1;
        while (true) {//将Xn,Yn装填到num[]中
            x = y / (a - (a - 1) * y);
            num[count][0] = x;

            if (num[count][0] < Xf) { //当Xn开始小于Xf时切换方程y为y‘
                y = a2 * x - b2;
            } else {//初始情况
                y = a1 * x + b1;
            }

            if (num[count][0] < Xw) { //当Xn开始小于Xw时停止循环
                break;
            }

            num[0][count + 1] = y;
            count++;
        }

        Xn = new double[20];
        Yn = new double[20];
        for (int i = 1; i <= count; i++) {
            Xn[i] = num[i][0];
            Yn[i] = num[0][i];
        }

        for (int i = 1; i < Xn.length; i++) {
            if (Xn[i] < Xf) {
                tbs = i - 1;
                djk = i;
                break;
            }
        }
        for (int i = 1; i < Xn.length; i++) {
            if (Xn[i] < Xw) {
                llb = i;
                break;
            }
        }

        System.out.println("精馏段有" + tbs + "块塔板，第" + djk + "块为进料板，全塔共有" + llb + "块理论板");

        for (int i = 1; i < Xn.length; i++) {
            if (Xn[i] == 0.0) {
                break;
            }
            System.out.print("\tY" + i + " = " + Yn[i]);
            System.out.println("\tX" + i + " = " + Xn[i] +
                    (Xn[i] < Xf ? " <Xf " : "")
                    +
                    (Xn[i] < Xw ? " <Xw " : "")
            );
        }

        HashMap<String, Double> map = new HashMap<>();
        map.put("yP", yP);
        map.put("Rm", Rm);
        map.put("R", R);
        map.put("L", L);
        map.put("V", V);
        map.put("L1", L1);
        map.put("V1", V1);
        map.put("a1", a1);
        map.put("b1", b1);
        map.put("a2", a2);
        map.put("b2", b2);
        map.put("llb", (double) llb);
        map.put("tbs", (double) tbs);
        return map;
    }


    Map<String, Double> chapter_3_2(double tbs, double llb) {
        System.out.println("**********3.2的计算开始**********");
        double Er = 0.52;
        System.out.println("全塔效率Er为0.52");
        int Nj, Nt;//精（提）馏段实际板数
        Nj = (int) (tbs / Er + 1);
        Nt = (int) ((llb - tbs) / Er + 1);
        HashMap<String, Double> map = new HashMap<>();
        map.put("Nj", (double) Nj);
        map.put("Nt", (double) Nt);
        return map;

    }


    Map<String, Double> chapter_4_1(double Nj, double Nt) {
        System.out.println("**********4.1的计算开始**********");
        //塔顶操作压力，每层塔板压降，进料板压力，精馏段平均压力，塔底操作压力，提馏段平均压力
        double PD, AP, PF, PM1, PW1, PM2;

        PD = P0 + Pg;
        AP = 0.7;
        PF = PD + Nj * AP;
        PM1 = (PD + PF) / 2;
        PW1 = PF + Nt * AP;
        PM2 = (PF + PW1) / 2;


        HashMap<String, Double> map = new HashMap<>();
        map.put("PD", PD);
        map.put("Ap", AP);
        map.put("PF", PF);
        map.put("PM1", PM1);
        map.put("PW1", PW1);
        map.put("PM2", PM2);
        return map;
    }


    Map<String, Double> chapter_4_2(double t1, double t2, double x1, double y1, double x2, double y2, double PM1, double PM2) {
        System.out.println("**********4.2的计算开始**********");
        double ML1, MV1, ML11, MV11;

        //pA,pB是t = t1时zww和zjw的密度 | pA1,pB1是t = t2时zww和zjw的密度
        double pA, pB, pA1, pB1;

        //pvm1精馏段气相平均密度 | pvm2提馏段气相平均密度
        double pvm1, pvm2;

        //pL1精馏段液相平均密度 | pV1提馏段液相平均密度
        double pL1, pV1;

        ML1 = zww * x1 + zjw * (1 - y1);
        MV1 = zww * y1 + zjw * (1 - x1);
        ML11 = zww * x2 + zjw * (1 - y2);
        MV11 = zww * y2 + zjw * (1 - x2);

        pA = quick_calculate(20,583.7,605.5,t1 - 40);
        pB = quick_calculate(20,620.0,638.9,t1 - 40);
        pA1 = quick_calculate(20,583.7,605.5,t2 - 40);
        pB1 = quick_calculate(20,620.0,638.9,t2 - 40);

        pvm1 = (PM1 * MV1) / (8.314 * (t1 + K));
        pvm2 = (PM2 * MV11) / (8.314 * (t2 + K));


        double tmp1 = (x1 * zww) / (x1 * zww + (1 - x1) * zjw);
        double tmp2 = (x2 * zww) / (x2 * zww + (1 - x2) * zjw);
        pL1 = 1 / (tmp1 / pA + (1 - tmp1) / pB);
        pV1 = 1 / (tmp2 / pA1 + (1 - tmp2) / pB1);

        HashMap<String, Double> map = new HashMap<>();
        map.put("ML1", ML1);
        map.put("MV1", MV1);
        map.put("ML11", ML11);
        map.put("MV11", MV11);
        map.put("pA", pA);
        map.put("pB", pB);
        map.put("pA1", pA1);
        map.put("pB1", pB1);
        map.put("pvm1", pvm1);
        map.put("pvm2", pvm2);
        map.put("pL1", pL1);
        map.put("pV1", pV1);
        return map;
    }


    Map<String, Double> chapter_4_3(double t1, double t2, double x1, double x2) {
        System.out.println("**********4.3的计算开始**********");
        double qa, qb, qm, qa1, qb1, qm1;

        qa = quick_calculate(20, 11.76, 13.85, t1 - 40);
        qb = quick_calculate(20, 13.23, 15.99, t1 - 40);
        qm = (qa * qb) / (qa * (1 - x1) + qb * x1);

        qa1 = quick_calculate(20, 11.76, 13.85, t2 - 40);
        qb1 = quick_calculate(20, 13.23, 15.99, t2 - 40);
        qm1 = (qa1 * qb1) / (qa1 * (1 - x2) + qb1 * x2);

        HashMap<String, Double> map = new HashMap<>();
        map.put("qa", qa);
        map.put("qb", qb);
        map.put("qm", qm);
        map.put("qa1", qa1);
        map.put("qb1", qb1);
        map.put("qm1", qm1);
        return map;

    }

    Map<String, Double> chapter_4_4(double x1, double x2) {
        System.out.println("**********4.4的计算开始**********");
        double Ua, Ub, Ua1, Ub1, U1, U2;

        Ua = quick_calculate(25, 7.37, 6.80, 15.485);
        Ub = quick_calculate(25, 7.10, 6.54, 15.485);
        Ua1 = quick_calculate(25, 7.96, 7.37, 5.965);
        Ub1 = quick_calculate(25, 7.66, 7.10, 6.3);
        U1 = x1 * Ua + Ub * (1 - x1);
        U2 = x2 * Ub1 + Ua1 * (1 - x2);

        HashMap<String, Double> map = new HashMap<>();
        map.put("Ua", Ua);
        map.put("Ub", Ub);
        map.put("Ua1", Ua1);
        map.put("Ub1", Ub1);
        map.put("U1", U1);
        map.put("U2", U2);
        return map;


    }


    void chapter_4_5() {
        System.out.println("**********4.5的计算开始**********");
        System.out.println("相对挥发度α="+a);
    }


    Map<String, Double> chapter_5_1(double pvm1, double pvm2, double pL1, double pV1, double L, double V, double ML1, double MV1, double ML11, double MV11, double l1, double v1) {
        System.out.println("**********5.1的计算开始**********");
        double PV1, PV2, PL1, PL2;

        PV1 = pvm1;
        PV2 = pvm2;
        PL1 = pL1;
        PL2 = pV1;

        double L1, V1, Ls1, Vs1, L2, V2, Ls2, Vs2;

        L1 = ML1 * L / 3600;
        V1 = MV1 * V / 3600;
        Ls1 = L1 / PL1;
        Vs1 = V1 / PV1;

        L2 = ML11 * l1 / 3600;
        V2 = MV11 * v1 / 3600;

        Ls2 = L2 / PL2;
        Vs2 = V2 / PV2;

        HashMap<String, Double> map = new HashMap<>();
        map.put("PV1", PV1);
        map.put("PV2", PV2);
        map.put("PL1", PL1);
        map.put("PL2", PL2);
        map.put("L1", L1);
        map.put("V1", V1);
        map.put("Ls1", Ls1);
        map.put("Vs1", Vs1);
        map.put("L2", L2);
        map.put("V2", V2);
        map.put("Ls2", Ls2);
        map.put("Vs2", Vs2);
        return map;
    }

    Map<String, Double> chapter_5_2(double Ls1, double Ls2, double Vs1, double Vs2, double PV1, double PV2, double PL1, double PL2, double qm, double qm1) {
        System.out.println("**********5.2的计算开始**********");

        System.out.println("精馏段>>>\n");

        //横坐标数值
        double hzbsz_1 = Ls1 / Vs1 * Math.pow((PL1 / PV1), 0.5);
        System.out.println("\t横坐标值为：" + hzbsz_1);

        double HT = 450;
        double hL = 60;
        //板间距
        double bjj = (HT - hL) / 1000;
        System.out.println("\tHT - hL = " + bjj);

        double C20_1 = 0.08;
        System.out.println("\t查表知C20=" + C20_1);

        double C = C20_1 * Math.pow((qm / 20), 0.2);
        System.out.println("\tC=" + C);

        double Vmax = C * Math.pow((PL1 - PV1) / PV1, 0.5);
        System.out.println("\tVmax=" + Vmax);

        //安全系数取0.8
        double aqxs = 0.8;
        System.out.println("\t安全系数取" + aqxs);

        double v1 = aqxs * Vmax;
        System.out.println("\tv1=" + v1);

        double D1 = Math.pow((4 * Vs1) / (3.14 * v1), 0.5);
        System.out.println("\tD1=" + D1);

        int tmp_1 = (int) (D1 * 10 + 1);
        double D1_after_1 = tmp_1 / 10.0;
        System.out.println("\tD1取整，D1=" + D1_after_1);

        double AT = 3.14 / 4 * D1 * D1;
        System.out.println("\tAT=" + AT);

        double v11 = Vs1 / AT;
        System.out.println("\tv1‘=" + v11);

        System.out.println("\n提馏段>>>\n");

        double hzbsz_2 = Ls2 / Vs2 * Math.pow((PL2 / PV2), 0.5);
        System.out.println("\t横坐标值为：" + hzbsz_2);

        //HT,hL,bjj与前面一样
        System.out.println("\tHT - hL = " + bjj);

        double C20_2 = 0.08;
        System.out.println("\t查表知C20=" + C20_2);

        double C1 = C20_2 * Math.pow((qm1 / 20), 0.5);
        System.out.println("\tC'=" + C1);

        double Vmax1 = C1 * Math.pow((PL2 - PV2) / PV2, 0.5);
        System.out.println("\tVmax'=" + Vmax);

        //aqxs不变
        System.out.println("\t安全系数取" + aqxs);

        double v2 = aqxs * Vmax1;
        System.out.println("\tv2=" + v2);

        double D2 = Math.pow((4 * Vs2) / (3.14 * v2), 0.5);
        System.out.println("\tD2=" + D2);

        int tmp_2 = (int) (D2 * 10 + 1);
        double D1_after_2 = tmp_2 / 10.0;
        System.out.println("\tD2取整，D2=" + D1_after_2);

        double AT1 = 3.14 / 4 * D2 * D2;
        System.out.println("\tAT'=" + AT1);

        double v21 = Vs2 / AT1;
        System.out.println("\tv2'=" + v21);


        HashMap<String, Double> map = new HashMap<>();
        map.put("D", D2);//取D1,D2大的数
        map.put("AT", AT1);
        map.put("hL", hL / 1000.0);
        map.put("HT", HT / 1000.0);
        return map;
    }


    Map<String, Double> chapter_5_3(double D, double AT, double Ls1, double hL, double HT, double Ls2) {
        System.out.println("**********5.3的计算开始**********");

        double lw, how, lA, hw, how1, hw1, Af, Wd, O, O1, v0, v01, h0, h01;

        System.out.println("(1) 堰长lw>>>");
        lw = 0.65 * D;
        System.out.println("\t取lw=" + lw);
        lA = Ls1;
        System.out.println("\t精馏段：");
        how = 2.84 / 1000 * Math.pow(lA * 3600 / lw, 2.0 / 3.0);
        hw = hL - how;
        System.out.println("\thow=" + how);
        System.out.println("\thw=" + hw);
        System.out.println("\t提馏段：");
        how1 = 2.84 / 1000 * Math.pow(Ls2 * 3600 / lw, 2.0 / 3.0);
        hw1 = hL - how1;
        System.out.println("\thow'=" + how1);
        System.out.println("\thw'=" + hw1);

        System.out.println("\n(2) 弓型降液管的宽度和横截面积>>>");
        System.out.println("\t查图得：Af / AT = 0.07\t Wd / D = 0.145");
        Af = 0.07 * AT;
        Wd = 0.145 * D;
        System.out.println("\t则：Af=" + Af + "\tWd=" + Wd);
        O = Af * HT / Ls1;
        System.out.println("\t精馏段：θ=" + O + (O > 5 ? " >5 " : " <5 "));
        O1 = Af * HT / Ls2;
        System.out.println("\t提馏段：θ'=" + O1 + (O1 > 5 ? " >5 " : " <5 "));

        System.out.println("\n(3) 降低管底隙高度>>>");
        System.out.println("v0=v0'=0.13");
        v0 = v01 = 0.13;
        h0 = Ls1 / (lw * v0);
        h01 = Ls2 / (lw * v01);
        System.out.println("\t精馏段：v0=" + v0 + "，h0=" + h0);
        System.out.println("\t提馏段：v0'=" + v01 + "，h0'=" + h01 + " (m) = " + h01 * 1000 + " (mm)");
        System.out.println("\t因为h0'不小20mm,故h0满足要求");


        HashMap<String, Double> map = new HashMap<>();
        map.put("Wd", Wd);
        map.put("Af", Af);
        map.put("hw", hw);
        map.put("hw1", hw1);
        map.put("lw", lw);
        map.put("h0", h0);
        map.put("h01", h01);
        return map;
    }

    Map<String, Double> chapter_5_4(double PV1, double PV2, double D, double Wd, double Vs1, double Vs2) {


        System.out.println("**********5.4的计算开始**********");
        System.out.println("(1) 塔板分布>>>");
        System.out.println("\t精馏段：[U0]Kp1=" + Math.pow(72.8 / PV1, 0.548));
        System.out.println("\t提馏段：[U0]Kp2=" + Math.pow(72.8 / PV2, 0.548));
        System.out.println("\t上下两段相应的阀孔动能因子为:");
        double F01 = Math.pow(72.8 / PV1, 0.548) * Math.pow(PV1, 0.5);
        System.out.println("\tF01=" + F01);
        System.out.println("\tF02=" + Math.pow(72.8 / PV2, 0.548) * Math.pow(PV2, 0.5));
        System.out.println("\t均属正常操作范围.");
        System.out.println("\n(2) 浮阀数目与排列>>>");
        double F0 = ((int)F01 + 1);
        System.out.println("\tF0="+ F0);
        System.out.println("\tv0=" + F0 / Math.pow(PV1, 0.5));
        System.out.println("\t精馏段：");
        System.out.println("\tWc=0.055m,Ws=0.065m");
        double R = D / 2 - 0.055;
        System.out.println("\tR=" + R);
        double x = D / 2 - Wd - 0.065;
        System.out.println("\tx=" + x);
        double Aa = 2 * (x * Math.pow(R * R - x * x, 0.5) + 3.14 * R * R / 180 * Math.toDegrees(Math.sin(x / R)));
        System.out.println("\tAa=" + Aa);
        System.out.println("\t提馏段：");
        System.out.println("\tWc=0.030m,Ws=0.055m");
        double R1 = D / 2 - 0.030;
        System.out.println("\tR=" + R1);
        double x1 = D / 2 - Wd - 0.055;
        System.out.println("\tx=" + x1);
        double Aa1 = 2 * (x1 * Math.pow(R1 * R1 - x1 * x1, 0.5) + (3.14 * R1 * R1 / 180 * Math.toDegrees(Math.asin(x1 / R1))));
        System.out.println("\tAa=" + Aa1);
        System.out.println("\n(3) 浮阀数n与开孔率>>>");

        double n_after = extracted_5_4("精馏段", F0, PV1, Vs1)[0];
        double n1_after = extracted_5_4("提馏段", F0, PV2, Vs2)[0];

        double n = extracted_5_4("精馏段", F0, PV1, Vs1)[1];
        double n1 = extracted_5_4("提馏段", F0, PV2, Vs2)[1];


        double t = 0.075;
        double ta = Aa / (n * t);
        double tb = Aa1 / (n1 * t);
        double tc = (int) ((Math.min(ta, tb)) * 100) / 100.0;
        System.out.println("\t取孔心距t=" + t);
        System.out.println("\t精馏段：t'=" + ta);
        System.out.println("\t提馏段：t'=" + tb);
        System.out.println("\t故取t'=" + tc + "m");
        System.out.println("\t重新计算孔速与阀数:");
        extracted_5_4_002("精馏段:", PV1, D, Vs1, Aa, n, t, ta);
        extracted_5_4_002("提馏段:", PV2, D, Vs2, Aa1, n, t, tb);

        HashMap<String, Double> map = new HashMap<>();
        map.put("N", n_after);
        map.put("N1", n1_after);
        return map;

    }

    void chapter_6_1(double D, double Wd, double AT, double Af, double PV1, double PV2, double PL1, double PL2) {
        System.out.println("**********6.1的计算开始**********");
        System.out.println("按泛点率80%计算");
        double Zl;
        double Ab;
        double a, b, c;
        double a1;
        System.out.println("\t精馏段：");
        Zl = D - 2 * Wd;
        System.out.println("\tZl=" + Zl);
        System.out.println("\t查物系数：K=1.0,Cf=0.125");
        Ab = AT - 2 * Af;
        System.out.println("\tAb=" + Ab);

        a = Math.pow(PV1 / (PL1 - PV1), 0.5);
        b = 1.36 * Zl;
        c = 0.8 * 0.125 * Ab;

        System.out.println("\t" + a + "Vs+" + b + "Ls=" + c);
        System.out.println("\tVs=" + (c / a) + "-" + (b / c) + "Ls");

        a1 = Math.pow(PV2 / (PL2 - PV2), 0.5);

        System.out.println("\t" + a1 + "Vs+" + b + "Ls=" + c);
        System.out.println("\tVs=" + (c / a1) + "-" + (b / c) + "Ls");
        System.out.println("精馏段\tLs(m3/s)\t0.002\t\t\t\t0.01");
        System.out.println("     \tVs(m3/s)\t" + ((c / a) - (b / c) * 0.002) + "\t" + ((c / a) - (b / c) * 0.01));
        System.out.println("提馏段\tLs(m3/s)\t0.002\t\t\t\t0.01");
        System.out.println("     \tVs(m3/s)\t" + ((c / a1) - (b / c) * 0.002) + "\t" + ((c / a1) - (b / c) * 0.01));
    }

    void chapter_6_2(double HT, double hw, double hw1, double PV1, double PV2, double N, double N1, double lw, double h0, double h01, double PL1, double PL2) {
        System.out.println("**********6.2的计算开始**********");
        extracted_6_2("精馏段", HT, hw, PV1, N, lw, h0, PL1);
        extracted_6_2("提馏段", HT, hw1, PV2, N1, lw, h01, PL2);
    }

    void chapter_6_3(double Af, double HT) {
        System.out.println("**********6.3的计算开始**********");
        System.out.println("(Ls)max=" + (Af * HT / 5));
    }

    void chapter_6_4(double N, double N1, double PV1, double PV2) {
        System.out.println("**********6.4的计算开始**********");
        System.out.println("(Vs1)min=" + Math.PI * Math.pow(0.039, 2) * N * 5 / Math.pow(PV1, 0.5));
        System.out.println("(Vs2)min=" + Math.PI * Math.pow(0.039, 2) * N1 * 5 / Math.pow(PV2, 0.5));
    }

    void chapter_6_5(double lw) {
        System.out.println("**********6.5的计算开始**********");
        double Lsmin = (0.006 * 1000 / 2.84) / Math.pow(3600 / lw, 2 / 3.0);
        double last = Math.pow(Math.pow(Lsmin, 3), 0.5);
        System.out.println("(Ls)min=" + last);
        System.out.println("**********所以计算均结束啦！！！**********");
    }

    private double[] extracted_5_4(String period,double F0,double PV,double Vs){
        System.out.println("\t"+period+"：");
        double u0 = F0 / Math.pow(PV, 0.5);
        System.out.println("\tu0=" + u0);
        double n = 4 * Vs / (u0 * 0.039 * 0.039 * 3.14);
        double n_after = (int) n + 1;
        System.out.println("\tn=" + n_after + "块");
        double fie = ((int) n + 1) * Math.pow(0.039, 2) / Math.pow(1.6, 2);
        System.out.println("\tφ=" + fie * 100 + "%");
        double[] doubles = new double[2];
        doubles[0] = n_after;
        doubles[1] = n;
        return doubles;
    }

    private void extracted_6_2(String period, double HT, double hw, double PV1, double N, double lw, double h0, double PL1) {
        double f;
        double e;
        double c;
        double d;
        double g;
        double b;
        double a;
        a = 0.5 * (HT + hw) - 1.5 * hw;
        b = 5.34 * PV1 / (Math.pow(Math.PI / 4 * 0.039 * 0.039 * N, 2) * 2 * 9.8 * PL1);
        c = 0.153 * 1 / Math.pow(lw * h0, 2);
        d = (1 + 0.5) * ((2.84 / 1000) * Math.pow(3600 / lw, 2 / 3.0));

        e = c / b * -1;
        f = d / b * -1;
        g = a / b;
        System.out.println("\t" + period + "：\tVs2=" + e + "Ls2 " + f + "Ls2/3 " + "+ " + g);

        System.out.println("\tLs1(m3/s)\t0.001\t\t\t\t0.003\t\t\t\t0.004\t\t\t\t0.007");
        System.out.println("\tVs1(m3/s)\t" +
                (e * Math.pow(0.001, 2) + f * Math.pow(0.001, 2 / 3.0) + g)
                + "\t" + (e * Math.pow(0.003, 2) + f * Math.pow(0.003, 2 / 3.0) + g)
                + "\t" + (e * Math.pow(0.004, 2) + f * Math.pow(0.004, 2 / 3.0) + g)
                + "\t" + (e * Math.pow(0.007, 2) + f * Math.pow(0.007, 2 / 3.0) + g)
        );


    }

    private void extracted_5_4_002(String period, double PV1, double D, double Vs1, double Aa, double n, double t, double ta) {
        System.out.println("\t" + period);
        int n_after_jl = (int) (Aa / (t * ta)) + 1;
        System.out.println("\tn=" + n_after_jl + "块");
        double u0_after_jl = 4 * Vs1 / (n * Math.pow(0.039, 2) * 3.14);
        System.out.println("\tu0=" + u0_after_jl);
        double F0_after_jl = u0_after_jl * Math.pow(PV1, 0.5);
        System.out.println("\tF0=" + F0_after_jl);
        double fie_after_jl = n_after_jl * Math.pow(0.039, 2) * Math.pow(D, 2);
        System.out.println("\tφ=" + fie_after_jl * 100 + "%");
    }

    /**
     * 计算相对挥发度
     */
    void calculate_relative_volatility(double t1, double t2) {
        double PA, PB, PA1, PB1, a1, a2;

        int up = 0,down;

        for (int i = 0; i < temp.length; i++) {
            if (t1 < temp[i]){
                up = i;
                break;
            }
        }
        down = up - 1;


        PA = quick_calculate(temp[up] - temp[down],array1[up],array1[down],t1 - temp[down]);
        PB = quick_calculate(temp[up] - temp[down],array2[up],array2[down],t1 - temp[down]);


        for (int i = 0; i < temp.length; i++) {
            if (t2 < temp[i]){
                up = i;
                break;
            }
        }
        down = up - 1;

        PA1 = quick_calculate(temp[up] - temp[down],array1[up],array1[down],t2 - temp[down]);
        PB1 = quick_calculate(temp[up] - temp[down],array2[up],array2[down],t2 - temp[down]);


        a1 = PA / PB;
        a2 = PA1 / PB1;

        a = Math.sqrt(a1 * a2);

        System.out.println("相对挥发度α=" + a);

    }

    /**
     * 快速计算
     * @param a 不做解释
     * @return 结果
     */
    double quick_calculate(double a, double b, double c, double d) {
        return d * (b - c) / a + c;
    }

    /**
     * 计算产品的摩尔分率
     */
    void calculate_mole_frac() {
        Xf = (mass_frac / zww) / (mass_frac / zww + (1.0 - mass_frac) / zjw);
        Xd = (product_concentration / zww) / (product_concentration / zww + (1.0 - product_concentration) / zjw);
        Xw = ((1 - product_concentration) / zww) / ((1 - product_concentration) / zww + (product_concentration) / zjw);
    }

    /**
     * 计算产品的平均摩尔质量
     *
     * @param mole_frac 摩尔分率
     * @return 平均摩尔质量
     */
    double calculate_average_molar_mass(double mole_frac) {
        return mole_frac * zww + (1 - mole_frac) * zjw;
    }

    /**
     * 物料衡算
     *
     * @param Mf 原料液的平均摩尔质量
     * @param Xw 塔釜产品的摩尔分率
     * @param Xf 原料液的摩尔分率
     * @param Xd 塔顶产品的摩尔分率
     */
    void material_balance(double Mf, double Xw, double Xf, double Xd) {
        F = year_process_mount * 1000 / (work_day * work_hours * Mf);

        W = (Xf - Xd) * F / (Xw - Xd);

        D = F - W;

    }

    /**
     * 初始化
     */
    void init(){
        System.out.println("全局变量如下：");
        System.out.println("年处理量：" + year_process_mount + "吨,进料组成：" + mass_frac * 100 + "%");

        temp = new double[]{36.1,40,45,50,55,60,65,68.7};
        array1 = new double[]{101.33,115.62,136.05,159.16,185.18,214.35,246.89,273.28};
        array2 = new double[]{31.98,37.26,45.02,54.05,64.66,76.36,89.96,101.33};

    }

    public static void main(String[] args) {

        HuaGong huaGong = new HuaGong();

        huaGong.init();

        huaGong.chapter_2_1();

        Map<String, Double> map = huaGong.chapter_2_2();
        System.out.println(map);
        Double t1 = map.get("t1");
        Double t2 = map.get("t2");
        Double x1 = map.get("x1");
        Double x2 = map.get("x2");
        Double y1 = map.get("y1");
        Double y2 = map.get("y2");

        Map<String, Double> map1 = huaGong.chapter_3_1(t1, t2);
        System.out.println(map1);
        Double L = map1.get("L");
        Double V = map1.get("V");
        Double L1 = map1.get("L1");
        Double V1 = map1.get("V1");

        Map<String, Double> map2 = huaGong.chapter_3_2(map1.get("tbs"), map1.get("llb"));
        System.out.println(map2);

        Map<String, Double> map3 = huaGong.chapter_4_1(map2.get("Nj"), map2.get("Nt"));
        System.out.println(map3);
        Double PM1 = map3.get("PM1");
        Double PM2 = map3.get("PM2");

        Map<String, Double> map4 = huaGong.chapter_4_2(t1, t2, x1, y1, x2, y2, PM1, PM2);
        System.out.println(map4);
        Double pvm1 = map4.get("pvm1");
        Double pvm2 = map4.get("pvm2");
        Double pL1 = map4.get("pL1");
        Double pV1 = map4.get("pV1");
        Double ML1 = map4.get("ML1");
        Double ML11 = map4.get("ML11");
        Double MV1 = map4.get("MV1");
        Double MV11 = map4.get("MV11");

        Map<String, Double> map5 = huaGong.chapter_4_3(t1, t2, x1, x2);
        System.out.println(map5);
        Double qm = map5.get("qm");
        Double qm1 = map5.get("qm1");

        Map<String, Double> map6 = huaGong.chapter_4_4(x1, x2);
        System.out.println(map6);

        huaGong.chapter_4_5();

        Map<String, Double> map7 = huaGong.chapter_5_1(pvm1, pvm2, pL1, pV1, L, V, ML1, MV1, ML11, MV11, L1, V1);
        System.out.println(map7);
        Double Ls1 = map7.get("Ls1");
        Double Ls2 = map7.get("Ls2");
        Double Vs1 = map7.get("Vs1");
        Double Vs2 = map7.get("Vs2");
        Double PL1 = map7.get("PL1");
        Double PV1 = map7.get("PV1");
        Double PV2 = map7.get("PV2");
        Double PL2 = map7.get("PL2");

        Map<String, Double> map8 = huaGong.chapter_5_2(Ls1, Ls2, Vs1, Vs2, PV1, PV2, PL1, PL2, qm, qm1);
        Double AT = map8.get("AT");
        Double D = map8.get("D");
        Double hL = map8.get("hL");
        Double HT = map8.get("HT");

        Map<String, Double> map9 = huaGong.chapter_5_3(D, AT, Ls1, hL, HT, Ls2);
        Double Wd = map9.get("Wd");
        Double Af = map9.get("Af");
        Double hw = map9.get("hw");
        Double hw1 = map9.get("hw1");
        Double lw = map9.get("lw");
        Double h0 = map9.get("h0");
        Double h01 = map9.get("h01");

        Map<String, Double> map10 = huaGong.chapter_5_4(PV1, PV2, D, Wd, Vs1, Vs2);
        Double N = map10.get("N");
        Double N1 = map10.get("N1");


        huaGong.chapter_6_1(D, Wd, AT, Af, PV1, PV2, PL1, PL2);

        huaGong.chapter_6_2(HT, hw, hw1, PV1, PV2, N, N1, lw, h0, h01, PL1, PL2);

        huaGong.chapter_6_3(Af, HT);

        huaGong.chapter_6_4(N, N1, PV1, PV2);

        huaGong.chapter_6_5(lw);

        HashMap<String, Double> map_total = new HashMap<>();
        System.out.println("打印所有的参数>>>");
        map_total.putAll(map);
        map_total.putAll(map1);
        map_total.putAll(map2);
        map_total.putAll(map3);
        map_total.putAll(map5);
        map_total.putAll(map6);
        map_total.putAll(map7);
        map_total.putAll(map8);
        map_total.putAll(map9);
        map_total.putAll(map10);
        System.out.println(map_total);
    }


}
