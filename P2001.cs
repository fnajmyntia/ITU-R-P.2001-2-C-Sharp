using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace AMS_IMT
{

    public class P2001
    {
        private string vDataFileSubdirectory;
        public string DataFileSubdirectory
        {
            get
            { return vDataFileSubdirectory; }
            set
            { vDataFileSubdirectory = value; }
        }

        private double[,] vTropoClim, vSurfwv_50_fixed, vEsarain_Pr6_v5, vEsarain_Mt_v5, vEsarain_Beta_v5;
        private double[,] vh0text, vDN_Median, vDN_SupSlope, vDN_SubSlope, vdndz_01, vFoEs50;
        private double[,] vFoEs10, vFoEs01, vFoEs0dot1;

        public double[,] TropoClim
        {
            get
            {
                if (vTropoClim == null)
                    vTropoClim = csvdouble2Dread("TropoClim.csv");
                return vTropoClim;
            }
        }

        public double[,] Surfwv_50_fixed
        {
            get
            {
                if (vSurfwv_50_fixed == null)
                    vSurfwv_50_fixed = csvdouble2Dread("Surfwv_50_fixed.csv");
                return vSurfwv_50_fixed;
            }
        }

        public double[,] Esarain_Pr6_v5
        {
            get
            {
                if (vEsarain_Pr6_v5 == null)
                    vEsarain_Pr6_v5 = csvdouble2Dread("Esarain_Pr6_v5.csv");
                return vEsarain_Pr6_v5;
            }
        }

        public double[,] Esarain_Mt_v5
        {
            get
            {
                if (vEsarain_Mt_v5 == null)
                    vEsarain_Mt_v5 = csvdouble2Dread("Esarain_Mt_v5.csv");
                return vEsarain_Mt_v5;
            }
        }

        public double[,] Esarain_Beta_v5
        {
            get
            {
                if (vEsarain_Beta_v5 == null)
                    vEsarain_Beta_v5 = csvdouble2Dread("Esarain_Beta_v5.csv");
                return vEsarain_Beta_v5;
            }
        }

        public double[,] h0text
        {
            get
            {
                if (vh0text == null)
                    vh0text = csvdouble2Dread("h0.csv");
                return vh0text;
            }
        }

        public double[,] DN_Median
        {
            get
            {
                if (vDN_Median == null)
                    vDN_Median = csvdouble2Dread("DN_Median.csv");
                return vDN_Median;
            }
        }

        public double[,] DN_SupSlope
        {
            get
            {
                if (vDN_SupSlope == null)
                    vDN_SupSlope = csvdouble2Dread("DN_SupSlope.csv");
                return vDN_SupSlope;
            }
        }

        public double[,] DN_SubSlope
        {
            get
            {
                if (vDN_SubSlope == null)
                    vDN_SubSlope = csvdouble2Dread("DN_SubSlope.csv");
                return vDN_SubSlope;
            }
        }

        public double[,] dndz_01
        {
            get
            {
                if (vdndz_01 == null)
                    vdndz_01 = csvdouble2Dread("dndz_01.csv");
                return vdndz_01;
            }
        }

        public double[,] FoEs50
        {
            get
            {
                if (vFoEs50 == null)
                    vFoEs50 = csvdouble2Dread("FoEs50.csv");
                return vFoEs50;
            }
        }

        public double[,] FoEs10
        {
            get
            {
                if (vFoEs10 == null)
                    vFoEs10 = csvdouble2Dread("FoEs10.csv");
                return vFoEs10;
            }
        }

        public double[,] FoEs01
        {
            get
            {
                if (vFoEs01 == null)
                    vFoEs01 = csvdouble2Dread("FoEs01.csv");
                return vFoEs01;
            }
        }

        public double[,] FoEs0dot1
        {
            get
            {
                if (vFoEs0dot1 == null)
                    vFoEs0dot1 = csvdouble2Dread("FoEs0.1.csv");
                return vFoEs0dot1;
            }
        }


        public P2001()
        {
            vDataFileSubdirectory = "";
        }

        public P2001(string P2001_DataFileDir)
        {
            vDataFileSubdirectory = P2001_DataFileDir;
        }

        private static readonly double[] probability_C21 = { 0.000555, 0.000802,    0.001139,    0.001594,    0.002196,
                0.002978,    0.003976,    0.005227,    0.006764,    0.008617,    0.010808,    0.013346,    0.016225,
                0.019419,    0.022881,    0.026542,    0.030312,    0.034081,    0.037724,    0.041110,    0.044104,
                0.046583,    0.048439,    0.049589,    0.049978,    0.049589,    0.048439,    0.046583,    0.044104,
                0.041110,    0.037724,    0.034081,    0.030312,    0.026542,    0.022881,    0.019419,    0.016225,
                0.013346,    0.010808,    0.008617,    0.006764,    0.005227,    0.003976,    0.002978,    0.002196,
                0.001594,    0.001139,    0.000802,    0.000555 };

        /// <summary>
        /// ITU-R P.2001-2,2015 RF Propagation Model
        /// A General Purpose Wide-Range Terrestrial Propagation Model in the Frequestion rang 30 MHz to 50 Ghz.
        /// Distance: 3-1,000km
        /// Antenna Altitudes: Up to 8000m above sea level
        /// </summary>
        /// <param name="freq">frequency (GHz): 30 MHz - 50 GHz</param>
        /// <param name="Tpol">Antenna Polarization: Vertical (1), Horizontal (0)</param>
        /// <param name="txlat">Latitude of Transmitter (degrees)</param>
        /// <param name="txlon">Longitude of Transmitter (degrees)</param>
        /// <param name="rxlat">Latitude of Receiver (degrees)</param>
        /// <param name="rxlon">Longitude of Receiver (degrees)</param>
        /// <param name="Hrg">Height of Electrical Center of Receiving Antenna above ground (meters)</param>
        /// <param name="Htg">Height of Electrical Center of Transmitting Antenna above ground (meters)</param>
        /// <param name="Tpc">Percentage of average year for which the prdicted basic transmission loss is not exceeded (%): Example 0.00001% to 99.99999%</param>
        /// <param name="Gr">Gain (dBi) of receiving antenna in the azimuthal direction of the path towards the other antenna, and at the elevation angle about the local horizontal of the other antenna in the case of a line-of-sight path, otherwise of the antenna's radio horizon, for median effective Earth radius.</param>
        /// <param name="Gt">Gain (dBi) of transmitting antenna in the azimuthal direction of the path towards the other antenna, and at the elevation angle about the local horizontal of the other antenna in the case of a line-of-sight path, otherwise of the antenna's radio horizon, for median effective Earth radius.</param>
        /// <param name="h">Elevation data (evenly spaced) in meters.</param>
        /// <returns></returns>
        public double p2001(double freq, double Tpol, double txlat, double txlon, double rxlat, double rxlon,
            double Hrg, double Htg, double Tpc, double Gr, double Gt, double[] h)
        {
            double Lb = 0.0;
            List<string> debug_table = new List<string>(); // This will be an out variable used for debugging
            const double lightspeed = 299792458.0; //Speed of Propagation: 299,792,458 m/s
            const double Re = 6371.0; //Average Earth Radius (km)
            const double erland = 22.0;  // Relative Permittivity for Land
            const double ersea = 80.0; // Relative Permittivity fo Sea
            const double sigmaland = 0.003; // Conducitivty for land (S/m)
            const double sigmasea = 5.0; // Conductivity for sea (S/m)
            //const bool SAVEDBUGDATA = false;  // Save Debug data

            double p_time, q_time, dn, hm;
            bool FlagShort;
            int mid_npts, mid_npts2, npts, mid_npts5;
            int profile_pts;
            double phi1qe, phi1qn, phi3qe, phi3qn, fsea;
            double dtm, dlm, dct, dcr, omega, h1, hh, hmid;
            double midlat, midlon;
            double[] t_lat, t_lon, tn_lat, tn_lon, d, d_tmp;
            double[] tropo_value;
            Vincenty vin = new Vincenty();
            // Reformat Data from .txt to .csv
            // Start of Section 3
            // Section 3.1 Limited Percentage Time
            p_time = Tpc + 0.00001 * ((50.0 - Tpc) / 50.0); //Equation 3.1.1
            q_time = 100.0 - p_time; //Equation 3.1.2
            //Great Circle Calculations: Attachment H or Matlab
            dn = vin.GetNewDist_m_SphericalEarth(txlat, txlon, rxlat, rxlon) / 1000.0; //Distance between tx and rx
            FlagShort = dn < .1; // If the distance is less than 100m or 0.1km.
            mid_npts = 3; // This is based upon the number of points in the terrain model (h).
            mid_npts5 = 5; // This is based upon the number of points in the terrain model (h).
            GenerateCoords(rxlat, rxlon, txlat, txlon, mid_npts, out t_lat, out t_lon, out d_tmp);  // PULL 3 Points
            midlat = t_lat[1]; midlon = t_lon[1];
            GenerateCoords(rxlat, rxlon, txlat, txlon, mid_npts5, out t_lat, out t_lon, out d_tmp);  // Pull 5 points

            phi3qe = t_lon[1];  // This is pulling the mid point
            phi3qn = t_lat[1];

            phi1qe = t_lon[3];  // This is pulling the 3/4 point
            phi1qn = t_lat[3];
            npts = h.Length; // This is based upon the number of points in the terrain model (h).
            GenerateCoords(rxlat, rxlon, txlat, txlon, npts, out t_lat, out t_lon, out d);   // Pull the Lat/Lon for each point  in h
            mid_npts2 = npts; // This is based upon the number of points in the terrain model (h).
            tn_lat = new double[t_lat.Length];
            tn_lon = new double[t_lon.Length];
            Array.Copy(t_lat, tn_lat, t_lat.Length);
            Array.Copy(t_lon, tn_lon, t_lon.Length);
            tropo_value = new double[tn_lat.Length];
            for (int i = 0; i < tn_lat.Length; ++i)
                tropo_value[i] = get_txt_value(tn_lat[i], tn_lon[i], TropoClim); // Gets the Climate of each point
            // At this point, calculate dtm/dml on tropo_value map, only worrying about
            // if it over land (1-6) or the sea (0) and using elevation data.
            find_ITU_dist(rxlat, rxlon, txlat, txlon, h, dn,
                out dtm, out dlm, out dct, out dcr, out omega);
            h1 = h[0];
            hh = h[h.Length - 1];
            profile_pts = h.Length;
            fsea = omega;
            hmid = h[(int)Math.Ceiling(h.Length / 2.0) - 1];
            double hts, hrs, hhi, hlo, ep, Sp, phirn, phire, phitn, phite;
            double phime, phimn, Nd65m1;
            // Start of Section 3.3 Antenna Altitude and Path Inclination
            hts = h[0] + Htg; //Equation 3.3.1a 
            hrs = h[h.Length - 1] + Hrg;  //Equation 3.3.1b

            //Assign the higher and lower antenna height above sea level
            hhi = Math.Max(hts, hrs); //Equation 3.3.2a
            hlo = Math.Min(hts, hrs); //Equation 3.3.2b

            //Calculate the Positive Value of Path Inclination:
            ep = (hhi - hlo) / dn; //Equation 3.3.3
            Sp = ep;
            // End of Section 3.3

            phitn = txlat;
            phite = txlon;
            phirn = rxlat;
            phire = rxlon;

            phime = midlon; // Longitude Point at Midpoint
            phimn = midlat; // Latitude Point at Midpoint

            Nd65m1 = get_txt_value(phimn, phime, dndz_01);

            double Sdn, Nd1km50, Nd1kmp, SdeltaNsup, SdeltaNsub, ae, Reff50, Reffp, ap, Cp;
            double thetae, lambda, dlt = 0.0;
            Sdn = get_txt_value(phimn, phime, DN_Median);
            Nd1km50 = -Sdn;
            SdeltaNsup = get_txt_value(phimn, phime, DN_SupSlope);
            SdeltaNsub = get_txt_value(phimn, phime, DN_SubSlope);

            if (p_time < 50)
                Nd1kmp = Nd1km50 + SdeltaNsup * Math.Log10(0.02 * p_time);  // Equation 3.4.1.2a
            else
                Nd1kmp = Nd1km50 - SdeltaNsub * Math.Log10(0.02 * q_time); //Equation 3.4.1.2b

            // End of Section 3.4
            // Start of Section 3.5 Effective Earth-radius geometry

            ae = (157 * Re) / (157 + Nd1km50); //Equation 3.5.1 (km), Median Effective Earth Radius
            Reff50 = ae;
            Cp = (157 + Nd1kmp) / (157 * Re);   //Equation 3.5.2 (km-1)
            //Effect Earth Curvature (cp is often positive, but it can be zero or negative.

            //Effective Earth Radius exceeded for p% time limited but to become infinite
            if (Cp > Math.Pow(10.0, -6.0))
                ap = 1 / Cp; //Equation 3.5.3a
            else
                ap = Math.Pow(10.0, 6.0); //Equation 3.5.3b
            Reffp = ap;
            // The path length expressed as the angle subtended by d km at the 
            // center of a sphere of effective Earth radius:
            thetae = dn / ae; // Original Equation 3.5.4
            // End of Section 3.5
            // Start of Section 3.6 Wavelength
            lambda = (Math.Pow(10.0, -9.0) * lightspeed) / freq; // Equation 3.6.1
            // End of Section 3.6
            // Start of Section 3.7: Path Classification and terminal horizon parameters
            // Find the highest elevation angle to an intermediate profile point, 
            // relative to the horzontal at the transmitter.
            double thetatim = double.MinValue, thetatr, thetat = 0.0, thetar = 0.0, dlr = 0.0, thetatpos, thetarpos;
            bool FlagLos50;
            int ithetatim = 0;
            for (int i = 2; i < d.Length - 2; ++i)
            {
                double tempval = ((h[i] - hts) / d[i]) - ((500.0 * d[i]) / ae);
                if (tempval > thetatim)   // Equation 3.7.1
                {
                    thetatim = tempval;
                    ithetatim = i;
                }
            }

            thetatr = ((hrs - hts) / dn) - ((500 * dn) / ae); // Equation 3.7.2

            // Case 1: Path is LOS 
            FlagLos50 = (thetatim < thetatr);
            int ilt = 0, ilr = 0;
            if (thetatim < thetatr) //The path is LOS
            {
                double vmax = double.MinValue, dim;
                int ivmax = 0;
                for (int i = 2; i < d.Length - 1; ++i)
                {
                    double tempvar = (h[i] + (((500 * d[i]) * (dn - d[i])) / ae) - ((hts * (dn - d[i]) + hrs * d[i]) / dn)) * Math.Sqrt((0.002 * dn) / (lambda * d[i] * (dn - d[i]))); //Equation 3.7.3 
                    if (tempvar > vmax)
                    {
                        tempvar = vmax;
                        ivmax = i;
                    }
                }
                //Find the itermediate profile point with the highest diffraction parameter
                dim = d[ivmax];  //distance where imax is the profile index which gives vmax
                dlt = dim;  //Equation 3.7.4a (km)
                dlr = dn - dim; //Equation 3.7.4b (km)
                ilt = ivmax; //3.7.4c  
                ilr = ivmax; //3.7.4d
                thetat = thetatr; //Equation 3.7.5a (mrad)
                thetar = -thetatr - ((1000 * dn) / ae); //Equation 3.7.5b (mrad)
            }
            // Case 2: Path is NLOS
            else if (thetatim >= thetatr) //The path is NLOS
            {
                int ithetarim = 0;
                double thetarim = double.MinValue;
                dlt = d[ithetatim]; //Equation 3.7.6a
                ilt = ithetatim; //Equation 3.7.6b
                thetat = thetatim; //Equation 3.7.7
                double[] temp_thetarim = new double[h.Length - 2];
                hm = 0.0;
                for (int i = 2; i < d.Length - 2; ++i)
                {
                    double tempvar = ((h[i] - hrs) / (dn - d[i])) - ((500 * (dn - d[i])) / ae);//Equation 3.7.8;
                    if (tempvar > thetarim)
                    {
                        thetarim = tempvar;
                        ithetarim = i; //i + 1;
                    }
                }
                dlr = dn - d[ithetarim];//Equation 3.7.9a
                ilr = ithetarim;//Equation 3.7.9b
                thetar = thetarim; //Equation 3.7.10
            }

            thetatpos = Math.Max(thetat, 0.0); // Equation 3.7.11a
            thetarpos = Math.Max(thetar, 0.0); // Equation 3.7.11b

            // End of Section 3.7
            // Start of Section 3.8 Effective Heights and Path Roughness Parameter
            double v1, v2, hstip, hsrip, hsripa, hstipa, mses, htea, hrea;
            v1 = 0.0;
            for (int i = 1; i < d.Length; ++i)
                v1 += (d[i] - d[i - 1]) * (h[i] + h[i - 1]); // Equation 3.8.1

            v2 = 0.0;
            for (int i = 1; i < d.Length; ++i)
                v2 += (d[i] - d[i - 1]) * (h[i] * (2 * d[i] + d[i - 1]) + h[i - 1] * (d[i] + 2 * d[i - 1])); //Equation 3.8.2

            hstip = (2 * v1 * dn - v2) / (dn * dn); //Equation 3.8.3a
            hsrip = (v2 - v1 * dn) / (dn * dn); //Equation 3.8.3b

            hstipa = Math.Min(hstip, h[0]); //Equation 3.8.4a
            hsripa = Math.Min(hsrip, h[h.Length - 1]); //3.8.4.b

            mses = (hsripa - hstipa) / dn; //Equation 3.8.5

            htea = hts - hstipa; //Equation 3.8.6a
            hrea = hrs - hsripa; //%Equation 3.8.6b;
            hm = double.MinValue;

            for (int i = ilt; i <= ilr; ++i)
                hm = Math.Max(hm, h[i] - (hstipa + mses * d[i])); // Equation 3.8.7

            double[] capH = new double[d.Length - 2];
            for (int i = 0; i < d.Length - 2; ++i)
                capH[i] = h[i + 1] - ((hts * (dn - d[i + 1]) + hrs * d[i + 1]) / dn); // Equation 3.8.8d

            double hobs = capH.Max(); //Equation 3.8.8a

            double aobr = double.MinValue, aobt = double.MinValue;
            double hst, hsr, gt, gr;

            for (int i = 0; i < d.Length - 2; ++i)
                aobt = Math.Max(aobt, capH[i] / d[i + 1]);

            for (int i = 0; i < d.Length - 2; ++i)
                aobr = Math.Max(aobr, capH[i] / (dn - d[i + 1]));

            if (hobs <= 0.0)
            {
                hst = hstip; // Equation 3.8.9a
                hsr = hsrip; // Equation 3.8.9b
            }
            else
            {
                gt = aobt / (aobt + aobr); // Equation 3.8.9e
                gr = aobr / (aobt + aobr); // Equation 3.8.9f
                hst = hstip - hobs * gt; // Equation 3.8.9c
                hsr = hsrip - hobs * gr; // Equation 3.8.9d
            }

            if (hst > h[0])
                hst = h[0]; // Equation 3.8.10a

            if (hsr > h[h.Length - 1])
                hsr = h[h.Length - 1]; // Equation 3.8.10b
            double htep, hrep, dtcv, drcv, hcv;
            htep = hts - hst; // Equation 3.8.11a
            hrep = hrs - hsr; // Equation 3.8.11b
            // End of Section 3.8
            // Start of Section 3.9 Tropospheric-scatter path segments
            dtcv = (dn * Math.Tan(0.001 * thetarpos + 0.5 * thetae) - 0.001 * (hts - hrs)) / (Math.Tan(0.001 * thetatpos + 0.5 * thetae) +
                Math.Tan(0.001 * thetarpos + 0.5 * thetae)); //Equation 3.9.1a
            if (dtcv > dn) // Limit dtcv to 0<=dtcv<=dn
                dtcv = dn;
            if (dtcv < 0) // Limit dtcv to 0<=dtcv<=dn
                dtcv = 0.0;
            drcv = dn - dtcv; // Equation 3.9.1b

            // Calculate phicve, phicvn by setting dpnt=dtcv in H.3.1 (Great Circle Path)(No need to use Appendix H)
            double Az, phicvn = 0.0, phicve = 0.0, phitcvn = 0.0, phitcve = 0.0, phircvn = 0.0, phircve = 0.0;
            double Aosur, Awsur, Awrsur, gamma0, gammaw, gammawr;

            Az = vin.GetBearing1to2_SphericalEarth(txlat, txlon, rxlat, rxlon);
            vin.GetNewSphericalEarthPoint(txlat, txlon, Az, dtcv, ref phicvn, ref phicve);

            hcv = hts + 1000 * dtcv * Math.Tan(0.001 * thetatpos) + ((1000 * dtcv * dtcv) / (2 * ae)); // Equation 3.9.2
            Az = vin.GetBearing1to2_SphericalEarth(txlat, txlon, rxlat, rxlon);
            vin.GetNewSphericalEarthPoint(txlat, txlon, Az, dtcv * 0.5, ref phitcvn, ref phitcve);

            // Calculate phitcve, phitcvn by setting dpnt=dn-0.5*drcv in H.3.1 (Great Circle Path)(No need to use Appendix H)
            vin.GetNewSphericalEarthPoint(txlat, txlon, Az, (dn - 0.5 * drcv), ref phircvn, ref phircve);
            vin.GetNewSphericalEarthPoint(txlat, txlon, Az, dtcv, ref phicvn, ref phicve);
            // End of Section 3.9
            // Start of Section 3.10 Gaseous Absorption on Surface Paths
            f2_att(freq, phimn, phime, hmid, hts, hrs, dn, out Aosur, out Awsur, out Awrsur, out gamma0, out gammaw, out gammawr);

            double Agsur = Aosur + Awsur; // Equation 3.10.1 The total gaseous attenuation under non-rain conditions

            // Section 3.11 Free-space basic transmission loss
            double LbfsD = 92.44 + 20 * Math.Log10(freq) + 20 * Math.Log10(dn); // Equation 3.11.1
            double Lbfs = LbfsD; // Equation 3.11.2
            // End of Section 3
            // Start of Section 4: Sub-models
            // There are 4 sub-models
            // Start of Section 4.1: Sub-model 1: Normal propagation close to the surface of the Earth
            // Use C.2 (Subroutine I)->B.2->B.4

            // Perform the preliminary rain/wet-snow calculation in C.2 with the following inputs
            double phie, phin, hrainlo, hrainhi, drain;
            phie = phime; // Equation 4.1.1a
            phin = phimn; // Equation 4.1.1b
            hrainlo = hlo; // Equation 4.1.1c
            hrainhi = hhi; // Equation 4.1.1d
            drain = dn;    //Equation 4.1.2

            double Fwvr, Q0ra, hRtop, Pr6, dr, kmod, alphamod, a1, b1, c1;
            double Q0ca, K, Q0cat, Q0car;
            double[] Pm, Gm;
            c2_calc_mergeC5(phin, phie, q_time, hrainlo, hrainhi, drain, Tpol, freq, hlo, hhi,
                out Fwvr, out Q0ra, out Pm, out Gm, out hRtop, out Pr6, out dr, out kmod, out alphamod, out a1, out b1, out c1); // Section C.2 Preliminary Calculation

            b2_Q0ca(Nd65m1, thetatim, thetatr, dn, ep, hlo, freq, phimn, dlt, thetat, hts, ilt, h, ilr, d, thetar, hrs, dlr,
                out Q0ca, out K, out Q0cat, out Q0car); // B.2 and B.3, Calculate Q0ca

            double A1 = I_Aiter_submodel1(Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1, q_time, Q0ca);

            double Ldbs, Ldsph, Ldbks, Ldba, Ldbka, Lbm1, Ld;
            double Lbm2, Lba, Aac, Aad, Aat;
            bool FlagLospa, FlagLosps;
            // Start of Section A: Diffraction Loss
            Ldsph = a2_loss(ap, htep, hrep, freq, dn, lambda, Tpol, omega); // Section A.2: Sphereical Earth Diffraction Loss
            a4_Ldba(d, h, Cp, hts, hrs, dn, lambda, out Ldba, out FlagLospa, out Ldbka); // Section A.4: Bullington Diffraction Loss for Actual Profile
            a5_Ldbs(d, dn, ap, htep, hrep, lambda, out Ldbs, out FlagLosps, out Ldbks); // A.5: Bullington Diffraction Loss for a notional smooth profile

            Ld = Ldba + Math.Max(Ldsph - Ldbs, 0.0); //Equation A.1.1 
            Lbm1 = Lbfs + Ld + A1 + Fwvr * (Awrsur - Awsur) + Agsur; //Equation 4.1.4 %OUTPUT of SUBMODEL 1
            // End of Section 4.1: Sub-model 1
            // Start of Section 4.2: Sub-model 2 Anomalous Propagation
            d_Lba(phimn, dlt, dlr, thetat, thetar, freq, ae, dn, hm, htea, hrea, p_time, q_time, dtm, dlm, dct, dcr, omega, hts, hrs,
                out Lba, out Aac, out Aad, out Aat); // Section D
            Lbm2 = Lba + Agsur; // Equation 4.2.1
            // End of Section 4.2
            // Start of Section 4.3: Sub-model 3:Troposcatter propagation
            double Lbs, thetas, Fwvrtx, Fwvrrx, A2t, A2r, A2;
            e_Lbs(phicvn, phicve, thetae, thetat, thetar, ae, dn, p_time, freq, Gt, Gr, Lbfs, out Lbs, out thetas); // Use Attachment E to calculate Lbs
            //Perform the preliminary rain/wet-snow calculation in C.2 for the TRANSMITTER to common-vloume path segment with the following inputs
            phie = phitcve; //Equation 4.3.1a
            phin = phitcvn; //Equation 4.3.1b
            hrainlo = hts; //Equation 4.3.1c
            hrainhi = hcv;    //Equation 4.3.1d
            drain = dtcv;    //Equation 4.3.1e

            // Save the value of Fwvr calculated in Section C.2 and call it Fwvrtx
            c2_calc_mergeC5(phin, phie, q_time, hrainlo, hrainhi, drain, Tpol, freq, hlo, hhi,
            out Fwvrtx, out Q0ra, out Pm, out Gm, out hRtop, out Pr6, out dr, out kmod, out alphamod, out a1, out b1, out c1); //Section C.2 Preliminary Calculation
            A2t = I_Aiter_submodel3(Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1, q_time);
            //Perform the preliminary rain/wet-snow calculation in C.2 for the RECEIVER to common-vloume path segment with the following inputs
            phie = phircve; //Equation 4.3.3a
            phin = phircvn; //Equation 4.3.3b
            hrainlo = hrs; //Equation 4.3.3c
            hrainhi = hcv;    //Equation 4.3.3d
            drain = drcv;    //Equation 4.3.3e
            // Save the value of Fwvr calculated in Section C.2 and call it Fwvrrx
            c2_calc_mergeC5(phin, phie, q_time, hrainlo, hrainhi, drain, Tpol, freq, hlo, hhi,
                out Fwvrrx, out Q0ra, out Pm, out Gm, out hRtop, out Pr6, out dr, out kmod, out alphamod, out a1, out b1, out c1); // Section C.2 Preliminary Calculation
            A2r = I_Aiter_submodel3(Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1, q_time);
            A2 = (A2t * (1 + 0.018 * dtcv) + A2r * (1 + 0.018 * drcv)) / (1 + 0.018 * dn); // Equation 4.3.6
            double Aos, Aws, Awrs, Aorcv, Aotcv, Awrcv, Awrrcv, Awrtcv, Awtcv, Ags, Lbm3, Lm, Lbm12;
            double Lbm4, Blank;
            f3_gas_absorption(h, phitn, phite, thetatpos, dtcv, phirn, phire, thetarpos, freq,
                out Aos, out Aws, out Awrs, out Aorcv, out Aotcv, out Awrcv, out Awrrcv, out Awrtcv, out Awtcv, ref drcv); // Use F.3 to Calculate Gaseous Attenuations Due to Oxygen
            Ags = Aos + Aws; //Equation 4.3.7
            Lbm3 = Lbs + A2 + 0.5 * (Fwvrtx + Fwvrrx) * (Awrs - Aws) + Ags; //Equation 4.3.8
            // End of Section 4.3
            // Start of Section 4.4: Sporadic-E
            double Lbe, foEs1, foEs2, Gamma1, Gamma2, LbEs1, LbEs2, Lp1r, Lp1t, Lp2r, Lp2t;
            g_sporadic(p_time, phimn, phime, dn, freq, ae, phi1qe, phi1qn, phi3qe, phi3qn, thetat, thetar, dlt, dlr,
                out Lbe, out foEs1, out foEs2, out Gamma1, out Gamma2, out LbEs1, out LbEs2, out Lp1r, out Lp1t, out Lp2r, out Lp2t); //Use Attachment G to calculate Lbe
            Lbm4 = Lbe; // Equation 4.4.1
            // End of Section 4.4
            // End of Section 4
            // Start of Section 5
            // Start of Section 5.1, Combining Sub-Models 1 and 2
            Lm = Math.Min(Lbm1, Lbm2);
            Lbm12 = Lm - 10 * Math.Log10(Math.Pow(10.0, -0.1 * (Lbm1 - Lm)) + Math.Pow(10.0, -0.1 * (Lbm2 - Lm))); // Equation 5.1.1
            // End of Section 5.1

            //Start of Section 5.2, Combining Sub Models 1+2, 3 and 4
            Lm = Math.Min(Math.Min(Lbm12, Lbm3), Lbm4);
            Lb = Lm - 5 * Math.Log10(Math.Pow(10.0, -0.2 * (Lbm12 - Lm)) + Math.Pow(10.0, -0.2 * (Lbm3 - Lm)) + Math.Pow(10.0, -0.2 * (Lbm4 - Lm))); // Equation 5.2.1
            // End of Section 5.2
            // End of Section 5

            // Save the debug data
            Blank = double.NaN;
            //if(SAVEDBUGDATA)
            //{
            //    string savefile = @"Debug.csv";
            //    StreamWriter sw = new StreamWriter(savefile);

            //    SaveDebug(sw, "FlagVP", FlagShort);
            //    SaveDebug(sw, "GHz", freq);
            //    SaveDebug(sw, "Grx", Gt);
            //    SaveDebug(sw, "Grt", Gr);
            //    SaveDebug(sw, "Hrg", Hrg);
            //    SaveDebug(sw, "Htg", Htg);
            //    SaveDebug(sw, "Phire", phire);
            //    SaveDebug(sw, "Phirn", phirn);
            //    SaveDebug(sw, "Phite", phite);
            //    SaveDebug(sw, "Phitn", phitn);
            //    SaveDebug(sw, "Tpc", Tpc);
            //    sw.WriteLine();
            //    SaveDebug(sw, "FlagLos50", FlagLos50);
            //    SaveDebug(sw, "FlagLospa", FlagLospa);
            //    SaveDebug(sw, "FlagLosps", FlagLosps);
            //    sw.WriteLine();
            //    sw.WriteLine();
            //    SaveDebug(sw, "A1", A1);
            //    SaveDebug(sw, "A2", A2);
            //    SaveDebug(sw, "A2r", A2r);
            //    SaveDebug(sw, "A2t", A2t);
            //    SaveDebug(sw, "Aac", Aac);
            //    SaveDebug(sw, "Aad", Aad);
            //    SaveDebug(sw, "Aat", Aat);
            //    SaveDebug(sw, "Ags", Ags);
            //    SaveDebug(sw, "Agsur", Agsur);
            //    SaveDebug(sw, "Aorcv", Aorcv);
            //    SaveDebug(sw, "Aos", Aos);
            //    SaveDebug(sw, "Aosur", Aosur);
            //    SaveDebug(sw, "Aotcv", Aotcv);
            //    SaveDebug(sw, "Awrcv", Awrcv);
            //    SaveDebug(sw, "Awrrcv", Awrrcv);
            //    SaveDebug(sw, "Awrs", Awrs);
            //    SaveDebug(sw, "Awrsur", Awrsur);
            //    SaveDebug(sw, "Awrtcv", Awrtcv);
            //    SaveDebug(sw, "Aws", Aws);
            //    SaveDebug(sw, "Awsur", Awsur);
            //    SaveDebug(sw, "Awtcv", Awtcv);
            //    sw.WriteLine();
            //    SaveDebug(sw, "Cp", Cp);
            //    SaveDebug(sw, "D", d);
            //    SaveDebug(sw, "Dcr", dcr);
            //    SaveDebug(sw, "Dct", dct);
            //    sw.WriteLine();
            //    SaveDebug(sw, "Dlm", dlm);
            //    SaveDebug(sw, "Dlr", dlr);
            //    SaveDebug(sw, "Dlt", dlt);
            //    SaveDebug(sw, "Drcv", drcv);
            //    SaveDebug(sw, "Dtcv", dtcv);
            //    SaveDebug(sw, "Dtm", dtm);
            //    SaveDebug(sw, "Foes1", foEs1);
            //    SaveDebug(sw, "Foes2", foEs2);
            //    SaveDebug(sw, "Fsea", fsea);
            //    SaveDebug(sw, "Fwvr", Fwvr);
            //    SaveDebug(sw, "Fwvrrx", Fwvrrx);
            //    SaveDebug(sw, "Fwvrtx", Fwvrtx);
            //    SaveDebug(sw, "GAM1", Gamma1);
            //    SaveDebug(sw, "GAM2", Gamma2);
            //    SaveDebug(sw, "Gamo", gamma0);
            //    SaveDebug(sw, "Gamw", gammaw);
            //    SaveDebug(sw, "Gamwr", gammawr);
            //    SaveDebug(sw, "H1", h1);
            //    SaveDebug(sw, "Hcv", hcv);
            //    SaveDebug(sw, "Hhi", hhi);
            //    SaveDebug(sw, "Hlo", hlo);
            //    SaveDebug(sw, "Hm", hm);
            //    SaveDebug(sw, "Hmid", hmid);
            //    SaveDebug(sw, "Hn?", 0);
            //    SaveDebug(sw, "Hrea", hrea);
            //    SaveDebug(sw, "Hrep", hrep);
            //    SaveDebug(sw, "Hrs", hrs);
            //    SaveDebug(sw, "Hsrip", hsrip);
            //    SaveDebug(sw, "Hsripa", hsripa);
            //    SaveDebug(sw, "Hstip", hstip);
            //    SaveDebug(sw, "Hstipa", hstipa);
            //    SaveDebug(sw, "Htea", htea);
            //    SaveDebug(sw, "Htep", htep);
            //    SaveDebug(sw, "Hts", hts);
            //    SaveDebug(sw, "Lb", Lb);
            //    SaveDebug(sw, "Lba", Lba);
            //    SaveDebug(sw, "Lbes1", LbEs1);
            //    SaveDebug(sw, "Lbes2", LbEs2);
            //    SaveDebug(sw, "Lbfs", Lbfs);
            //    SaveDebug(sw, "Lbm1", Lbm1);
            //    SaveDebug(sw, "Lbm2", Lbm2);
            //    SaveDebug(sw, "Lbm3", Lbm3);
            //    SaveDebug(sw, "Lbm4", Lbm4);
            //    SaveDebug(sw, "Lbs", Lbs);
            //    SaveDebug(sw, "Ld", Ld);
            //    SaveDebug(sw, "Ldba", Ldba);
            //    SaveDebug(sw, "Ldbka", Ldbka);
            //    SaveDebug(sw, "Ldbks", Ldbks);
            //    SaveDebug(sw, "Ldbs", Ldbs);
            //    SaveDebug(sw, "dLdsph", Ldsph);
            //    SaveDebug(sw, "Lp1r", Lp1r);
            //    SaveDebug(sw, "Lp1t", Lp1t);
            //    SaveDebug(sw, "Lp2r", Lp2r);
            //    SaveDebug(sw, "Lp2t", Lp2t);
            //    SaveDebug(sw, "Mses", mses);
            //    SaveDebug(sw, "N", profile_pts);
            //    SaveDebug(sw, "Nd1km50", Nd1km50);
            //    SaveDebug(sw, "Nd1kmp", Nd1kmp);
            //    SaveDebug(sw, "Nd65m1", Nd65m1);
            //    SaveDebug(sw, "Nlr", ilr);
            //    SaveDebug(sw, "Nlt", ilt);
            //    sw.WriteLine();
            //    sw.WriteLine();
            //    sw.WriteLine();
            //    sw.WriteLine();
            //    SaveDebug(sw, "Phi1qe", phi1qe);
            //    SaveDebug(sw, "Phi1qn", phi1qn);
            //    SaveDebug(sw, "Phi3qe", phi3qe);
            //    SaveDebug(sw, "Phi3qn", phi3qn);
            //    SaveDebug(sw, "Phicve", phicve);
            //    SaveDebug(sw, "Phicvn", phicvn);
            //    SaveDebug(sw, "Phime", phime);
            //    SaveDebug(sw, "Phimn", phimn);
            //    SaveDebug(sw, "Phircve", phircve);
            //    SaveDebug(sw, "Phircvn", phircvn);
            //    SaveDebug(sw, "Phitcve", phitcve);
            //    SaveDebug(sw, "Phitcvn", phitcvn);
            //    SaveDebug(sw, "Qoca", Q0ca);
            //    SaveDebug(sw, "Reff50", Reff50);
            //    SaveDebug(sw, "Reffp", Reffp);
            //    SaveDebug(sw, "Sp", Sp);
            //    SaveDebug(sw, "Thetae", thetae);
            //    SaveDebug(sw, "Thetar", thetar);
            //    SaveDebug(sw, "Thetarpos", thetarpos);
            //    SaveDebug(sw, "Thetas", thetas);
            //    SaveDebug(sw, "Thetat", thetat);
            //    SaveDebug(sw, "Thetatpos", thetatpos);
            //    SaveDebug(sw, "Tpcp", p_time);
            //    SaveDebug(sw, "Tpcq", q_time);
            //    //sw.WriteLine();
            //    //SaveDebug(sw, "Wave", wave);
            //    //SaveDebug(sw, "Wvsur", wvsur);
            //    //SaveDebug(sw, "WvSurrx", wvSurrx);
            //    //SaveDebug(sw, "WvSurtx", wvSurtx);
            //    //SaveDebug(sw, "Ztropo", ztropo);

            //    sw.Close();
            //}
            return Lb;
        }

        void SaveDebug(StreamWriter sw, string s, object v)
        {
            sw.Write(string.Format("{0},", s)); sw.WriteLine(v);
        }

        // P.2001-2-2015: Attachment G3: 2 Hop Propagation
        void g3_LbEs2(double p_time, double phi1qe, double phi1qn, double phi3qe, double phi3qn, double dn,
            double freq, double ae, double thetat, double thetar, double dlt, double dlr,
            out double LbEs2, out double foEs2h, out double Gamma2, out double Lp2r, out double Lp2t)
        {
            double foEs21, foEs22, hes, l2, Lbfs2, alpha2;
            double er2, d2t, d2r, v2t, v2r;
            foEs21 = g1_foEs(p_time, phi1qn, phi1qe); // Obtain foEs for the quarter point of the path
            foEs22 = g1_foEs(p_time, phi3qn, phi3qe); // Obtain foEs for the quarter point of the path
            foEs2h = Math.Min(foEs21, foEs22);
            Gamma2 = ((40.0 / (1.0 + (dn / 260.0) + Math.Pow(dn / 500.0, 2.0))) + 0.2 * Math.Pow(dn / 5200.0, 2.0)) * Math.Pow((1000.0 * freq) / foEs2h, 2.0) + Math.Exp((dn - 3220.0) / 560.0); // Equation G.3.1
            hes = 120; // hes is the height of the sporadic-E layer in km, set to 120km.
            l2 = 4.0 * Math.Pow(Math.Pow(ae, 2.0) + Math.Pow(ae + hes, 2.0) - 2 * ae * (ae + hes) * Math.Cos(dn / (4 * ae)), 0.5); // Equation G.3.2
            Lbfs2 = 92.44 + 20 * Math.Log10(freq) + 20.0 * Math.Log10(l2); // %Equation 3.11.1 %Equation G.3.3, free-space loss calculated for the slope distance.
            alpha2 = dn / (4 * ae); // Equation G.3.4a
            er2 = 0.5 * Math.PI - Math.Atan((ae * Math.Sin(alpha2)) / (hes + ae * (1 - Math.Cos(alpha2)))) - alpha2; // Equation G.3.4
            d2t = 0.001 * thetat - er2; // Equation G.3.5
            d2r = 0.001 * thetar - er2; // Equation G.3.5
            if (d2t >= 0.0)
                v2t = 3.651 * Math.Sqrt(1000 * freq * dlt * ((1 - Math.Cos(d2t)) / (Math.Cos(0.001 * thetat)))); // Equation G.3.6a
            else
                v2t = -3.651 * Math.Sqrt(1000 * freq * dlt * ((1 - Math.Cos(d2t)) / (Math.Cos(0.001 * thetat)))); // Equation G.3.6b 
            if (d2r >= 0.0)
                v2r = 3.651 * Math.Sqrt(1000 * freq * dlr * ((1 - Math.Cos(d2r)) / (Math.Cos(0.001 * thetar)))); // Equation G.3.6a
            else
                v2r = -3.651 * Math.Sqrt(1000 * freq * dlr * ((1 - Math.Cos(d2r)) / (Math.Cos(0.001 * thetar)))); // Equation G.3.6b 
            Lp2t = knife_edge_Jv(v2t); // Equation G.3.7a
            Lp2r = knife_edge_Jv(v2r); // Equation G.3.7b
            LbEs2 = Lbfs2 + Gamma2 + Lp2t + Lp2r; // Equation G.3.8
        }

        // P.2001-2-2015: Attachment G4: Basic Transmission Loss
        void g_sporadic(double p_time, double phimn, double phime, double dn, double freq, double ae, double phi1qe, double phi1qn,
            double phi3qe, double phi3qn, double thetat, double thetar, double dlt, double dlr,
            out double Lbe, out double foEs1, out double foEs2, out double Gamma1, out double Gamma2, out double LbEs1,
            out double LbEs2, out double Lp1r, out double Lp1t, out double Lp2r, out double Lp2t)
        {
            // This method should not be considered reliable at low or high
            // geomagnetic latitudes, and it need not be calculated for a LoS path.

            // It should be noted that incidents of high signal strength due to this
            // phenomenon exhibit a very strong seasonal dependence.
            g2_LbEs1(p_time, phimn, phime, dn, freq, ae, thetat, thetar, dlt, dlr, out LbEs1, out foEs1, out Gamma1, out Lp1r, out Lp1t);
            g3_LbEs2(p_time, phi1qe, phi1qn, phi3qe, phi3qn, dn, freq, ae, thetat, thetar, dlt, dlr, out LbEs2, out foEs2, out Gamma2, out Lp2r, out Lp2t);
            Lbe = -10 * Math.Log10(Math.Pow(10.0, -0.1 * LbEs1) + Math.Pow(10.0, -0.1 * LbEs2)); // Equation G.4.1c
            if (LbEs1 < (LbEs2 - 20.0))
            {
                Lbe = LbEs1; // Equation G.4.1a
            }
            if (LbEs2 < (LbEs1 - 20.0))
            {
                Lbe = LbEs2; // Equation G.4.1b
            }
        }

        // P.2001-2-2015: Attachment G2: 1 Hop Propagation
        void g2_LbEs1(double p_time, double phimn, double phime, double dn, double freq, double ae,
            double thetat, double thetar, double dlt, double dlr,
            out double LbEs1, out double foEs1, out double Gamma1, out double Lp1r, out double Lp1t)
        {
            double foEs, hes, l1, Lbfs1, alpha1, er1;
            double d1t, d1r, v1t, v1r;
            foEs = g1_foEs(p_time, phimn, phime); // Obtain foEs for the midpoint of the path.
            foEs1 = foEs;
            Gamma1 = ((40.0 / (1 + (dn / 130.0) + Math.Pow(dn / 250, 2.0))) + 0.2 * Math.Pow(dn / 2600, 2.0)) * Math.Pow((1000 * freq) / foEs, 2.0) + Math.Exp((dn - 1660) / 280); // Equation G.2.1
            hes = 120; // hes is the height of the sporadic-E layer in km, set to 120km.
            l1 = 2 * Math.Pow(Math.Pow(ae, 2.0) + Math.Pow(ae + hes, 2.0) - 2 * ae * (ae + hes) * Math.Cos(dn / (2 * ae)), 0.5); // Equation G.2.2
            Lbfs1 = 92.44 + 20 * Math.Log10(freq) + 20.0 * Math.Log10(l1); // %Equation 3.11.1 // Equation G.2.3, free-space loss calculated for the slope distance.
            alpha1 = dn / (2.0 * ae); // Equation G.2.4a
            er1 = 0.5 * Math.PI - Math.Atan((ae * Math.Sin(alpha1)) / (hes + ae * (1.0 - Math.Cos(alpha1)))) - alpha1; // Equation G.2.4
            d1t = 0.001 * thetat - er1; // Equation G.2.5
            d1r = 0.001 * thetar - er1; // Equation G.2.5
            if (d1t >= 0.0)
                v1t = 3.651 * Math.Sqrt(1000 * freq * dlt * ((1 - Math.Cos(d1t)) / (Math.Cos(0.001 * thetat)))); // Equation G.2.6a
            else
                v1t = -3.651 * Math.Sqrt(1000 * freq * dlt * ((1 - Math.Cos(d1t)) / (Math.Cos(0.001 * thetat)))); // Equation G.2.6b 
            if (d1r >= 0.0)
                v1r = 3.651 * Math.Sqrt(1000 * freq * dlr * ((1 - Math.Cos(d1r)) / (Math.Cos(0.001 * thetar)))); // Equation G.2.6a
            else
                v1r = -3.651 * Math.Sqrt(1000 * freq * dlr * ((1 - Math.Cos(d1r)) / (Math.Cos(0.001 * thetar)))); // Equation G.2.6b 
            Lp1t = knife_edge_Jv(v1t); // Equation G.2.7a
            Lp1r = knife_edge_Jv(v1r); // Equation G.2.7b
            LbEs1 = Lbfs1 + Gamma1 + Lp1t + Lp1r; // Equation G.2.8
        }

        // P.2001-2-2015: Attachment G1: Derivation of foEs
        double g1_foEs(double p_time, double lat, double lon)
        {
            // This method should not be considered reliable at low or high
            // geomagnetic latitudes, and it need not be calculated for a LoS path.

            // It should be noted that incidents of high signal strength due to this
            // phenomenon exhibit a very strong seasonal dependence.
            double foEs, p1 = 0.0, p2 = 0.0, foEs1 = 0.0, foEs2 = 0.0;
            // Table G.1
            if (p_time < 1.0)
            {
                p1 = 0.1;
                p2 = 1;
                foEs1 = get_txt_value(lat, lon, FoEs0dot1);
                foEs2 = get_txt_value(lat, lon, FoEs01);
            }
            else if (p_time >= 1.0 && p_time <= 10.0)
            {
                p1 = 1.0;
                p2 = 10.0;
                foEs1 = get_txt_value(lat, lon, FoEs01);
                foEs2 = get_txt_value(lat, lon, FoEs10);
            }
            else if (p_time > 10.0)
            {
                p1 = 10;
                p2 = 50;
                foEs1 = get_txt_value(lat, lon, FoEs10);
                foEs2 = get_txt_value(lat, lon, FoEs50);
            }
            foEs = foEs1 + (foEs2 - foEs1) * ((Math.Log10(p_time / p1)) / Math.Log10(p2 / p1)); // Equation G.1.1 (MHz)
            return foEs;
        }

        //P.2001-2-2015: Attachment F.3 Gaseous Absorption for a Troposcatter Path
        void f3_gas_absorption(double[] h, double phitn, double phite, double thetatpos, double dtcv, double phirn, double phire, double thetarpos, double freq,
            out double Aos, out double Aws, out double Awrs, out double Aorcv, out double Aotcv,
            out double Awrcv, out double Awrrcv, out double Awrtcv, out double Awtcv, ref double drcv)
        {
            double hsur, dcv, thetaelev;
            hsur = h[0];
            thetaelev = thetatpos;
            dcv = dtcv;
            f4_path(hsur, thetaelev, dcv, phitn, phite, freq, out Aotcv, out Awtcv, out Awrtcv); //Use Equation F.3.1a-c
            hsur = h.Length;
            thetaelev = thetarpos;
            dcv = drcv;
            f4_path(hsur, thetaelev, dcv, phirn, phire, freq, out Aorcv, out Awrcv, out Awrrcv); //Use Equation F.3.2a-c
            Aos = Aotcv + Aorcv; //Equation F.3.3a
            Aws = Awtcv + Awrcv; //Equation F.3.3b
            Awrs = Awrtcv + Awrrcv; //Equation F.3.3c
        }

        //P.2001-2-2015: Attachment F.4 Gaseous Absorption for a Terminal/Common-Volume Troposcatter Path
        void f4_path(double hsur, double thetaelev, double dcv, double lat, double lon, double freq, out double Ao, out double Aw, out double Awr)
        {
            double dox, dw, deo, dew;
            double gamma0, gammaw, gammawr;
            f6_att(freq, lat, lon, hsur, out gamma0, out gammaw, out gammawr); //Use F.6.2
            dox = 5 / (0.65 * Math.Sin(0.001 * thetaelev) + 0.35 * Math.Sqrt((Math.Pow(Math.Sin(0.001 * thetaelev), 2.0)) + 0.00304)); //Equation F.4.1a
            dw = 2 / (0.65 * Math.Sin(0.001 * thetaelev) + 0.35 * Math.Sqrt((Math.Pow(Math.Sin(0.001 * thetaelev), 2.0)) + 0.00122)); //Equation F.4.1b
            deo = dox * (1 - Math.Exp(-1 * dcv / dox)) * Math.Exp(-1 * hsur / 5000); //Equation F.4.2a
            dew = dw * (1 - Math.Exp(-1 * dcv / dw)) * Math.Exp(-1 * hsur / 2000); //Equation F.4.2.b
            Ao = gamma0 * deo; //Equation F.4.3a
            Aw = gammaw * dew; //Equation F.4.3b
            Awr = gammawr * dew; //Equation F.4.3c
        }

        //P.2001-2-2015: Attachement I: Iterative Procedure to invert a cumulative distribution function.
        //It seems like there should be a way to do this with Matlab
        double I_Aiter_submodel3(double Q0ra, double Pr6, double[] Pm, double[] Gm, double hRtop,
            double hrainlo, double dr, double kmod, double alphamod, double a1, double b1, double c1, double q_time)
        {
            double Aiter, Ainit, Ahigh, Alow, Astep, qhigh, qlow, q;
            //Section I.2
            //Stage 1, setting the search range
            Ainit = 10.0; //10dB
            Ahigh = Ainit / 2; //Equation I.2.1
            Alow = -Ainit / 2; //Equation I.2.2
            Astep = Ainit; //Equation I.2.3
            qhigh = Qiter_sm3(Ahigh, Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1); //Equation I.2.4a
            qlow = Qiter_sm3(Alow, Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1); //Equation I.2.4b
            q = q_time;
            //Stage 1 
            for (int i = 1; i <= 10; ++i) //No more than 10 loops
            {
                if (q < qhigh)
                {
                    Alow = Ahigh;
                    qlow = qhigh;
                    Astep = 2 * Astep;
                    Ahigh = Ahigh + Astep;
                    qhigh = Qiter_sm3(Ahigh, Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1);
                }
                if (q > qlow)
                {
                    Ahigh = Alow;
                    qhigh = qlow;
                    Astep = 2 * Astep;
                    Alow = Alow - Astep;
                    qlow = Qiter_sm3(Alow, Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1);
                }
                if (q >= qhigh && q <= qlow)
                    break; //Proceed to Stage 2
            }

            //Stage 2: Binary Search
            double Atry, qtry, Aacc;
            int niter;
            Atry = 0.5 * (Alow + Ahigh); //Equation I.2.5
            qtry = Qiter_sm3(Atry, Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1);  //Equation I.2.6
            Aacc = 0.01;
            niter = (int)Math.Ceiling(3.32 * Math.Log10(Astep / Aacc));
            for (int i = 1; i <= niter; ++i)
            {
                if (qtry < q)
                    Ahigh = Atry;
                else
                    Alow = Atry;
                Atry = 0.5 * (Alow + Ahigh); //Equation I.2.5
                qtry = Qiter_sm3(Atry, Q0ra, Pr6, Pm, Gm, hRtop, hrainlo, dr, kmod, alphamod, a1, b1, c1);  //Equation I.2.6
            }
            Aiter = Atry;
            return Aiter;
        }

        //P.2001-2-2015: Attachment E: Troposcatter
        void e_Lbs(double phicvn, double phicve, double thetae, double thetat, double thetar, double ae,
            double dn, double p_time, double freq, double Gt, double Gr, double Lbfs,
            out double Lbs, out double thetas)
        {
            double mp_trop_cli, htrop, He, LN, ds, Y90 = 0.0, C, Yp;
            mp_trop_cli = get_txt_value(phicvn, phicve, TropoClim);

            double[] M = new double[] { 129.60, 119.73, 109.30, 128.50, 119.73, 123.2, 116.0 }; //From Table E.1
            double[] gamma = new double[] { 0.33, 0.27, 0.32, 0.27, 0.27, 0.27, 0.27 }; //From Table E.1
            thetas = 1000 * thetae + thetat + thetar; //Equation E.1
            htrop = 0.125 * Math.Pow(10.0, -6.0) * Math.Pow(thetas, 2.0) * ae; //Equation E.4
            He = 0.25 * Math.Pow(10, -3.0) * dn; //Equation E.3
            if (mp_trop_cli == 0.0)
                LN = 20 * Math.Log10(5 + gamma[6] * He) + 4.34 * gamma[6] * htrop; //Equation E.2
            else
                LN = 20 * Math.Log10(5 + gamma[(int)mp_trop_cli - 1] * He) + 4.34 * gamma[(int)mp_trop_cli - 1] * htrop; //Equation E.2
            ds = 0.001 * thetas * ae; //Equation E.5
            switch ((int)mp_trop_cli)
            {
                case 0:
                    Y90 = -9.5 - 3 * Math.Exp(-0.137 * htrop); //Equation E.7
                    break;
                case 1:
                    if (ds < 100.0)
                        Y90 = -8.2;//Equation E.8
                    else if (ds < 1000.0 && ds >= 100.0)
                        Y90 = 1.006 * Math.Pow(10.0, -8.0) * Math.Pow(ds, 3.0) - 2.569 * Math.Pow(10.0, -5.0) * Math.Pow(ds, 2.0) + 0.2242 * ds - 10.2; //Equation E.8
                    else if (ds >= 1000.0)
                        Y90 = -3.4; //Equation E.8
                    break;
                case 3:
                    if (ds < 100.0)
                        Y90 = -10.845; //Equation E.9
                    else if (ds < 465.0 && ds >= 100.0)
                        Y90 = -4.5 * Math.Pow(10.0, -7.0) * Math.Pow(ds, 3.0) + 4.45 * Math.Pow(10.0, -4.0) * Math.Pow(ds, 2.0) - 0.122 * ds - 2.645; //Equation E.9
                    else if (ds >= 465.0)
                        Y90 = -8.4; //Equation E.9
                    break;
                case 4:
                    if (ds < 100.0)
                        Y90 = -11.5;//Equation E.10
                    else if (ds < 550.0 && ds >= 100.0)
                        Y90 = -8.519 * Math.Pow(10.0, -8.0) * Math.Pow(ds, 3.0) + 7.444 * Math.Pow(10.0, -5.0) * Math.Pow(ds, 2.0) - 4.18 * Math.Pow(10.0, -4.0) * ds - 12.1; //Equation E.10
                    else if (ds >= 550.0)
                        Y90 = -4.0; //Equation E.10
                    break;
                case 2:
                case 5:
                case 6:
                    Y90 = -2.2 - (8.1 - 0.23 * Math.Min(freq, 4.0)) * Math.Exp(-0.137 * htrop); //Equation E.6
                    break;
            }
            if (p_time >= 50.0)
                C = 1.26 * Math.Pow(-1 * Math.Log10((100 - p_time) / 50), 0.63); //Euqation E.11a
            else
                C = -1.26 * Math.Pow(-1 * Math.Log10(p_time / 50), 0.63); //Equation E.11b
            Yp = C * Y90; //Equation E.12
            if (thetas < Math.Pow(10.0, -6.0)) //Limit the value of thetas to thetas>=10^-6
                thetas = Math.Pow(10.0, -6.0);
            double Ldist, Lfreq, Lcoup;
            Ldist = Math.Max(10.0 * Math.Log10(dn) + 30.0 * Math.Log10(thetas) + LN, 20.0 * Math.Log10(dn) + 0.573 * thetas + 20.0); //Equation E.13
            Lfreq = 25 * Math.Log10(freq) - 2.5 * Math.Pow(Math.Log10(0.5 * freq), 2.0); //Equation E.14
            Lcoup = 0.07 * Math.Exp(0.055 * (Gt + Gr)); //Equation E.15
            if (mp_trop_cli == 0.0)
                Lbs = M[6] + Lfreq + Ldist + Lcoup - Yp; //Equation E.16
            else
                Lbs = M[(int)mp_trop_cli - 1] + Lfreq + Ldist + Lcoup - Yp; //Equation E.16
            //To Avoid under-estimating troposcatter loss for short paths, limits Lbs
            if (Lbs < Lbfs)
                Lbs = Lbfs; //Equation E.17 
        }


        void d_Lba(double phimn, double dlt, double dlr, double thetat, double thetar, double freq, double ae, double dn,
            double hm, double htea, double hrea, double p_time, double q_time, double dtm, double dlm, double dct, double dcr, double omega,
            double hts, double hrs,
            out double Lba, out double Aac, out double Aad, out double Aat)
        {
            double tao, u1, u2, u3, u4, Beta0, gtr, grr, Alf, gammad;
            double thetast, thetasr, Ast, Asr, Act, Acr, Betaduct;
            double thetaat, thetaar, thetaa, dar, alpha, CapGamma;
            //Start of D.2: Point Incidence of ducting
            tao = Math.Round(1 - Math.Exp(-0.000412 * (Math.Pow(dlm, 2.41))), 6); //Equation D.2.1
            u1 = Math.Pow(Math.Pow(10.0, ((-1 * dtm) / (16 - 6.6 * tao))) + Math.Pow(10.0, (-1 * (2.48 + 1.77 * tao))), 0.2); //Equation D.2.2
            if (u1 > 1) //Limit the value of u1<=1
                u1 = 1;
            if (Math.Abs(phimn) <= 70.0)
                u4 = Math.Pow(10.0, ((-0.935 + 0.0176 * Math.Abs(phimn)) * Math.Log10(u1))); //Equation D.2.3
            else
                u4 = Math.Pow(10.0, (0.3 * Math.Log10(u1))); //Equation D.2.3
            if (Math.Abs(phimn) <= 70.0)
                Beta0 = (Math.Pow(10.0, (-0.015 * Math.Abs(phimn) + 1.67))) * u1 * u4; //Equation D.2.4
            else
                Beta0 = 4.17 * u1 * u4; //Equation D.2.4
            //End of D.2
            //Start of D.3: Site-Shielding losses with respect to anomalous propagation mechanisms
            gtr = 0.1 * dlt; //Equation D.3.1a
            grr = 0.1 * dlr; //Equation D.3.1b
            thetast = thetat - gtr; //Equation D.3.2a
            thetasr = thetar - grr; //Equation D.3.2b
            if (thetast > 0.0)
                Ast = 20.0 * Math.Log10(1 + 0.361 * thetast * Math.Pow(freq * dlt, 0.5)) + 0.264 * thetast * Math.Pow(freq, 1.0 / 3.0); //Equation D.3.3a
            else
                Ast = 0.0; //Equation D.3.3b

            if (thetasr > 0.0)
                Asr = 20.0 * Math.Log10(1 + 0.361 * thetasr * Math.Pow(freq * dlr, 0.5)) + 0.264 * thetasr * Math.Pow(freq, 1.0 / 3.0); //Equation D.3.4a
            else
                Asr = 0.0;  //Equation D.3.4b
            //End of D.3

            //Start of D.4: Over-sea surface duct coupling corrections
            if (omega >= 0.75 && dct <= dlt && dct <= 5.0) //5km
                Act = -3.0 * Math.Exp(-0.25 * Math.Pow(dct, 2.0)) * (1 + Math.Tanh(0.07 * (50.0 - hts))); //Equation D.4.2a
            else
                Act = 0.0; //Equation D.4.2b

            if (omega >= 0.75 && dcr <= dlr && dcr <= 5) //5km
                Acr = -3.0 * Math.Exp(-0.25 * Math.Pow(dcr, 2.0)) * (1 + Math.Tanh(0.07 * (50.0 - hrs))); //Equation D.4.3a
            else
                Acr = 0.0;  //Equation D.4.3b

            //End of D.4

            //Start of D.5
            if (freq < 0.5)
                Alf = (45.375 - 137.0 * freq + 92.5 * Math.Pow(freq, 2.0)) * omega; //Equation D.5.2a
            else
                Alf = 0.0;  //Equation D.5.2b
            Aac = 102.45 + 20.0 * Math.Log10(freq * (dlt + dlr)) + Alf + Ast + Asr + Act + Acr; //Equation D.5.1
            //End of D.5

            //Start of D.6
            gammad = 5.0 * Math.Pow(10.0, -5.0) * ae * (Math.Pow(freq, 1.0 / 3.0)); //Equation D.6.1
            thetaat = Math.Min(thetat, gtr); //Equation D.6.2a
            thetaar = Math.Min(thetar, grr); //Equation D.6.2b
            thetaa = ((1000.0 * dn) / ae) + thetaat + thetaar; //Equation D.6.3
            Aad = gammad * thetaa; //Equation D.6.4 (Correct EQ)
            //End of D.6

            //Start of D.7
            dar = Math.Min(dn - dlt - dlr, 40.0); //Equation D.7.1
            if (hm > 10.0)
                u3 = Math.Exp(-4.6 * Math.Pow(10.0, -5.0) * (hm - 10.0) * (43.0 + (6.0 * dar))); //Equation D.7.2a
            else
                u3 = 10.0; //Equation D.7.2b
            alpha = -0.6 - 3.5 * Math.Pow(10.0, -9.0) * Math.Pow(dn, 3.1) * tao; //Equation D.7.3

            if (alpha < -3.4)
                alpha = -3.4;
            u2 = Math.Pow((500 * Math.Pow(dn, 2.0)) / (ae * Math.Pow(Math.Sqrt(htea) + Math.Sqrt(hrea), 2.0)), alpha); // Equation D.7.4
            if (u2 > 1.0)
                u2 = 1.0;
            Betaduct = Beta0 * u2 * u3; //D.7.5
            CapGamma = (1.076 * Math.Exp(Math.Pow(-10.0, -6.0) * Math.Pow(dn, 1.13) * (9.51 - 4.8 * Math.Log10(Betaduct) + 0.198 * Math.Pow(Math.Log10(Betaduct), 2.0)))) /
                (Math.Pow(2.0058 - Math.Log10(Betaduct), 1.012)); //Equation D.7.6
            Aat = -12 + ((1.2 + 0.0037 * dn) * Math.Log10(p_time / Betaduct)) + (12.0 * Math.Pow(p_time / Betaduct, CapGamma)) + (50.0 / q_time); //Equation D.7.7
            //End of D.7
            Lba = Aac + Aad + Aat; // Equation D.8.1
        }

        // P.2001-2-2015: Attachement A.5: Bullington Diffraction Loss for a notional smooth profile
        void a5_Ldbs(double[] d, double dn, double ap, double htep, double hrep, double lambda, out double Ldbs, out bool FlagLosps, out double Ldbks)
        {
            double Stim = double.MinValue;
            for (int i = 1; i < d.Length - 1; ++i)
                Stim = Math.Max((500 * (dn - d[i]) / ap) - (htep / d[i]), Stim); //Equation A.5.1
            double vs, db, vb;
            double Str;
            Str = (hrep - htep) / dn; //Equation A.5.2 
            FlagLosps = (Stim < Str);
            Ldbks = 0.0;
            if (Stim < Str) //Case 1: LOS 
            {
                vs = double.MinValue;
                for (int i = 1; i < d.Length - 1; ++i)
                {
                    double vs1, vs2, vs3;
                    vs1 = (500 * d[i] * (dn - d[i])) / ap;
                    vs2 = (htep * (dn - d[i]) + (hrep * d[i])) / dn;
                    vs3 = Math.Sqrt((0.002 * dn) / (lambda * d[i] * (dn - d[i])));
                    vs = Math.Max((vs1 - vs2) * vs3, vs); // Equation A.5.3
                }
                db = double.NaN;
                vb = double.NaN;
                Ldbks = knife_edge_Jv(vs); //Equation A.5.4
            }
            else if (Stim >= Str) //Case 2: NLOS 
            {
                double Srim = double.MinValue;
                for (int i = 1; i < d.Length - 1; ++i)
                    Srim = Math.Max(((500 * d[i]) / ap) - (hrep / (dn - d[i])), Srim); //Equation A.5.5
                vs = double.NaN;
                db = (hrep - htep + Srim * dn) / (Stim + Srim); //Equation A.5.6
                vb = ((htep + Stim * db) - ((htep * (dn - db) + hrep * db) / dn)) * Math.Sqrt((0.002 * dn) / (lambda * db * (dn - db))); //Equation A.5.7
                Ldbks = knife_edge_Jv(vb); //Equation A.5.8
            }
            Ldbs = Ldbks + (1.0 - Math.Exp(-Ldbks / 6.0)) * (10.0 + 0.02 * dn); //Equation A.5.9
        }

        // P.2001-2-2015: Attachement A.4: Bullington Diffraction loss for actual profile
        void a4_Ldba(double[] d, double[] h, double cp, double hts, double hrs, double dn, double lambda,
            out double Ldba, out bool FlagLospa, out double Ldbka)
        {
            double vb, db;
            Ldbka = 0.0;
            double Stim = double.MinValue, Srim = double.MinValue, va = double.MinValue;
            for (int i = 1; i < d.Length - 2; ++i)
                Stim = Math.Max(Stim, (h[i] + (500.0 * cp * d[i] * (dn - d[i])) - hts) / d[i]); //Equation A.4.1 
            double Str;
            Str = (hrs - hts) / dn; // Equation A.4.2

            FlagLospa = (Stim < Str);

            if (Stim < Str) //Case 1: LOS 
            {
                for (int i = 1; i < d.Length - 2; ++i)
                    va = Math.Max(va, (h[i] + 500 * cp * d[i] * (dn - d[i]) - ((hts * (dn - d[i]) + hrs * d[i]) / dn)) * Math.Sqrt((0.002 * dn) / (lambda * d[i] * (dn - d[i])))); // Equation A.4.3  
                vb = double.NaN;
                db = double.NaN;
                Srim = double.NaN;
                Ldbka = knife_edge_Jv(va); // Equation A.4.4  
            }
            else if (Stim >= Str) // Case 2: NLOS 
            {
                for (int i = 1; i < d.Length - 2; ++i)
                    Srim = Math.Max((h[i] + 500 * cp * d[i] * (dn - d[i]) - hrs) / (dn - d[i]), Srim); // Equation A.4.5 {Correct}
                db = (hrs - hts + (Srim * dn)) / (Stim + Srim); // Equation A.4.6 {Correct}
                va = double.NaN;
                vb = (hts + (Stim * db) - ((hts * (dn - db) + (hrs * db)) / dn)) * Math.Sqrt((0.002 * dn) / (lambda * db * (dn - db))); // Equation A.4.7 {Correct}
                Ldbka = knife_edge_Jv(vb); // Equation A.4.8
            }
            Ldba = Ldbka + (1 - Math.Exp((-1 * Ldbka) / 6)) * (10 + (0.02 * dn)); // Equation A.4.9 {Correct}
        }

        // P.2001-2-2015: 3.12: Knife-Edge Diffraction Loss
        double knife_edge_Jv(double v)
        {
            double Jv;
            if (v > -0.78)
                Jv = 6.9 + (20.0 * Math.Log10(Math.Sqrt(Math.Pow(v - 0.1, 2.0) + 1.0) + v - 0.1)); // Equation 3.12.1a 
            else
                Jv = 0.0; // Equation 3.12.1b
            return Jv;
        }

        // P.2001-2-2015: Attachement I: Iterative Procedure to invert a cumulative distribution function.
        double I_Aiter_submodel1(double Q0ra, double Pr6, double[] Pm, double[] Gm, double hRtop, double hrainlo, double dr,
            double kmod, double alphamod, double a1, double b1, double c1, double q_time, double Q0ca)
        {
            double Aiter, Ainit, Ahigh, Alow, Astep, qhigh, qlow, q;
            // Section I.2
            // Stage 1, setting the search range
            Ainit = 10.0; // 10dB
            Ahigh = Ainit / 2.0; // Equation I.2.1
            Alow = -Ainit / 2.0; // Equation I.2.2
            Astep = Ainit; // Equation I.2.3
            qhigh = Qiter_sm1(Ahigh, Q0ca, Q0ra, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1); // Equation I.2.4a
            qlow = Qiter_sm1(Alow, Q0ca, Q0ra, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1); // Equation I.2.4b
            q = q_time;

            // Stage 1 
            for (int i = 1; i <= 10; ++i) // No more than 10 loops
            {
                if (q < qhigh)
                {
                    Alow = Ahigh;
                    qlow = qhigh;
                    Astep = 2.0 * Astep;
                    Ahigh = Ahigh + Astep;
                    qhigh = Qiter_sm1(Ahigh, Q0ca, Q0ra, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1);
                }
                if (q > qlow)
                {
                    Ahigh = Alow;
                    qhigh = qlow;
                    Astep = 2.0 * Astep;
                    Alow = Alow - Astep;
                    qlow = Qiter_sm1(Alow, Q0ca, Q0ra, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1);
                }
                if (q >= qhigh && q <= qlow)
                    break; // Proceed to Stage 2
            }

            // Stage 2: Binary Search
            double Atry, qtry, Aacc;
            int niter;
            Atry = 0.5 * (Alow + Ahigh); // Equation I.2.5
            qtry = Qiter_sm1(Atry, Q0ca, Q0ra, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1);  // Equation I.2.6
            Aacc = 0.01;
            niter = (int)Math.Ceiling(3.32 * Math.Log10(Astep / Aacc));
            for (int i = 1; i <= niter; ++i)
            {
                if (qtry < q)
                    Ahigh = Atry;
                else
                    Alow = Atry;
                Atry = 0.5 * (Alow + Ahigh); // Equation I.2.5
                qtry = Qiter_sm1(Atry, Q0ca, Q0ra, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1);  // Equation I.2.6
            }
            Aiter = Atry;
            return Aiter;
        }

        // P.2001-2-2015: Qiter(q) for Submodel 1
        double Qiter_sm1(double A, double Q0ca, double Q0ra, double hrainlo, double hRtop, double[] Pm, double[] Gm,
            double Pr6, double dr, double kmod, double alphamod, double a1, double b1, double c1)
        {
            double varout, Qcaf, Qrain;
            Qcaf = b4_Qcaf(Q0ca, A); // B.4: Calculate Qcaf(A)
            Qrain = c3_Qrain(A, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1); // C.3: Calculate Qrain(A)
            varout = Qrain * (Q0ra / 100) + Qcaf * (1 - (Q0ra / 100)); // Qiter(A):Equation 4.1.3
            return varout;
        }

        // P.2001-2-2015: Attachement B.4: Percentage time a given clear-air fade level is exceeded on a surface path
        double b4_Qcaf(double Q0ca, double A)
        {
            double Qcaf = 0.0, qt, qa, qs, qe;
            if (A >= 0.0)
            {
                qt = 3.576 - 1.955 * Math.Log10(Q0ca); // Equation B.4.1.b
                qa = 2.0 + (1.0 + 0.3 * Math.Pow(10.0, -0.05 * A)) * Math.Pow(10.0, -0.016 * A) * (qt + 4.3 * (Math.Pow(10.0, -0.05 * A) + (A / 800.0))); // Equation B.4.1a
                Qcaf = 100.0 * (1.0 - Math.Exp(-Math.Pow(10.0, -0.05 * qa * A) * Math.Log(2.0))); // B.4.1 %ln(2)
            }
            else if (A < 0.0)
            {
                qs = -4.05 - 2.35 * Math.Log10(Q0ca); // Equation B.4.2b
                qe = 8.0 + (1.0 + 0.3 * Math.Pow(10.0, 0.05 * A)) * (Math.Pow(10.0, 0.035 * A)) * (qs + 12.0 * (Math.Pow(10.0, 0.05 * A) - (A / 800.0))); // Equation B.4.2a
                Qcaf = 100.0 * Math.Exp(-Math.Pow(10.0, 0.05 * qe * A) * Math.Log(2.0)); // Equation B.4.2  %ln(2)
            }
            return Qcaf;
        }

        // P.2001-2-2015: Attachement B.2: Characterize Multi-Path Activity
        void b2_Q0ca(double Nd65m1, double thetatim, double thetatr, double dn, double ep, double hlo, double freq, double phimn, double dlt, double thetat, double hts, int ilt, double[] h, int ilr, double[] d, double thetar, double hrs, double dlr,
            out double Q0ca, out double K, out double Q0cat, out double Q0car)
        {
            double dca, eca, hca, dcat, ecat, hcat, dcar, ecar, hcar;
            Q0ca = 0.0;
            Q0cat = 0.0;
            Q0car = 0.0;
            K = Math.Pow(10.0, (-1 * (4.6 + 0.0027 * Nd65m1))); // Equation B.2.1
            // Case 1: Path is LOS 
            if (thetatim < thetatr) // The path is LOS
            {
                dca = dn; // Equation B.2.2a
                eca = ep; // Equation B.2.2.b
                hca = hlo; // Equation B.2.2.c
                Q0car = double.NaN;
                Q0cat = double.NaN;
                Q0ca = b3_Q0ca(dca, eca, hca, freq, phimn, K); // USE B.3
            }

            // Case 2: Path is NLOS
            else if (thetatim >= thetatr) // The path is NLOS
            {
                dcat = dlt; // Equation B.2.3a
                ecat = Math.Abs(thetat); // Equation B.2.3b
                hcat = Math.Min(hts, h[ilt]); // Equation B.2.3c
                Q0cat = b3_Q0ca(dcat, ecat, hcat, freq, phimn, K); // USE B.3 to Calculate Q0cat

                // Calculate Q0car
                dcar = dlr; // Equation B.2.4a
                // dcar=46.348
                ecar = Math.Abs(thetar); // Equation B.2.4b
                hcar = Math.Min(hrs, h[ilr]); // Equation 4.2c
                Q0car = b3_Q0ca(dcar, ecar, hcar, freq, phimn, K); // USE B.3 to Calculate Q0car
                Q0ca = Math.Max(Q0cat, Q0car); // Equation B.2.5
            }
        }

        void f2_att(double freq, double phimn, double phime, double hmid, double hts, double hrs, double dn,
            out double Aosur, out double Awsur, out double Awrsur, out double gamma0, out double gammaw, out double gammawr)
        {
            double hsur = hmid; // The terrain heigh at the middle of the path
            f6_att(freq, phimn, phime, hsur, out gamma0, out gammaw, out gammawr);
            double hrho = 0.5 * (hts + hrs); // Equation F.2.1
            Aosur = gamma0 * dn * Math.Exp(-hrho / 5000); // Equation F.2.2a
            Awsur = gammaw * dn * Math.Exp(-hrho / 2000); // Equation F.2.2b
            Awrsur = gammawr * dn * Math.Exp(-hrho / 2000); // Equation F.2.2c
        }

        // P.2001-2-2015: Attachment F6(F5): Gaseous Absoprtion for Surface Path
        // Not valid for frequencies greater than 54 GHz. P2001 is limited to 50 GHz
        // See ITU-R P.676 for a more general expression
        void f6_att(double freq, double lat, double lon, double hsur, out double gamma0, out double gammaw, out double gammawr)
        {
            // F.6 Specific Sea Level Attenuations
            gamma0 = ((7.2 / (Math.Pow(freq, 2.0) + 0.34)) + (0.62 / (Math.Pow(54.0 - freq, 1.16) + 0.83))) * Math.Pow(freq, 2.0) * Math.Pow(10, -3.0); // Equation F.6.1
            double psur = get_txt_value(lat, lon, Surfwv_50_fixed);
            double psea = psur * Math.Exp(hsur / 2000); // Equation F.6.2.b
            double n1 = 0.955 + 0.006 * psea; // Equation F.6.2a
            gammaw = (0.046 + 0.0019 * psea + (3.98 * n1) / (Math.Pow(freq - 22.235, 2.0) + 9.42 * Math.Pow(n1, 2.0)) *
                (1 + Math.Pow((freq - 22) / (freq + 22), 2.0))) * Math.Pow(freq, 2.0) * psea * Math.Pow(10, -4.0); // Equation F.6.2
            double psurr, psur2, psear, nr;
            if (hsur <= 2600.00)
                psurr = psur + 0.4 + 0.0003 * hsur; // Equation F.5.1
            else
                psurr = psur + 5 * Math.Exp(-hsur / 1800); // Equation F.5.1
            psur2 = psurr; // Re-evaluate psur
            psear = psur2 * Math.Exp(hsur / 2000); // Equation F.6.2.b
            nr = 0.955 + 0.006 * psear; // Equation F.6.2a
            gammawr = (0.046 + 0.0019 * psear + (3.98 * nr) / (Math.Pow(freq - 22.235, 2.0) + 9.42 * Math.Pow(nr, 2.0)) *
                (1 + Math.Pow((freq - 22) / (freq + 22), 2.0))) * Math.Pow(freq, 2.0) * psear * Math.Pow(10, -4.0); // Equation F.6.2
        }


        void GenerateCoords(double lat1, double lon1, double lat2, double lon2, int npts, out double[] lat, out double[] lon, out double[] d)
        {
            double dist_m, spacing_m, bearing, curdist_m = 0.0;
            int ptr = 0;
            Vincenty vin = new Vincenty();
            dist_m = vin.GetNewDist_m_SphericalEarth(lat1, lon1, lat2, lon2);
            bearing = vin.GetBearing1to2_SphericalEarth(lat1, lon1, lat2, lon2);
            spacing_m = dist_m / (double)(npts - 1);
            lat = new double[npts];
            lon = new double[npts];
            d = new double[npts];
            lat[0] = lat1;
            lon[0] = lon1;
            d[0] = 0.0;
            do
            {
                ++ptr;
                d[ptr] = spacing_m * (double)ptr;
                vin.GetNewSphericalEarthPoint(lat1, lon1, bearing, d[ptr], ref lat[ptr], ref lon[ptr]);
                d[ptr] /= 1000.0;  // convert to km
            } while (ptr < npts - 1);
        }

        public double[,] csvdouble2Dread(string fname)
        {
            string FullFileName = Path.Combine(DataFileSubdirectory, fname);
            double[,] data;
            List<double> lstdata = new List<double>();
            int r = 0, c;
            using (StreamReader sr = new StreamReader(FullFileName))
            {
                string s;
                while (sr.Peek() >= 0)
                {
                    s = sr.ReadLine();
                    string[] p = s.Split(',');
                    for (int i = 0; i < p.Length; ++i)
                        lstdata.Add(double.Parse(p[i]));
                    ++r;
                }
            }
            c = lstdata.Count() / r;
            data = new double[r, c];
            int ptr = 0;
            for (int rr = 0; rr < r; ++rr)
                for (int cc = 0; cc < c; ++cc)
                    data[rr, cc] = lstdata[ptr++];
            return data;
        }

        //P.2001-2-2015: Attachement A.2: Spherical-Earch Diffraction Loss
        double a2_loss(double ap, double htep, double hrep, double freq, double dn, double lambda, double Tpol, double omega)
        {
            double Ldsph, dlos, adft, Ldft, csph, msph, bsph;
            double d1, d2, hsph, aem;
            dlos = Math.Sqrt(2 * ap) * (Math.Sqrt(0.001 * htep) + Math.Sqrt(0.001 * hrep)); //Equation A.2.1
            if (dn >= dlos)
            {
                adft = ap;
                //Calculate diffraction Loss using A.3 to find Ldft
                Ldft = a3_loss(freq, adft, Tpol, htep, hrep, omega, dn);  //Subrountine for A.3
                Ldsph = Ldft;
            }
            else
            {
                csph = (htep - hrep) / (htep + hrep); //Equation A.2.2d
                msph = (250 * Math.Pow(dn, 2.0)) / (ap * (htep + hrep)); //Equation A.2.2e
                bsph = 2.0 * Math.Sqrt((msph + 1.0) / (3.0 * msph)) * Math.Cos(Math.PI / 3.0 + 1.0 / 3.0 * Math.Acos((3.0 * csph / 2.0) *
                    Math.Sqrt((3.0 * msph) / (Math.Pow(msph + 1, 3.0))))); //Equation A.2.2c
                d1 = dn / 2 * (1 + bsph); //Equation A.2.2.a
                d2 = dn - d1; //Equation A.2.2b
                hsph = ((htep - 500 * (Math.Pow(d1, 2.0) / ap)) * d2 + (hrep - 500 * (Math.Pow(d2, 2.0) / ap)) * d1) / dn; //Equation A.2.2
                hrep = 17.456 * Math.Sqrt((d1 * d2 * lambda) / dn); //Equation A.2.3
                if (hsph > hrep)
                    Ldsph = 0;
                else
                {
                    aem = 500 * Math.Pow((dn / (Math.Sqrt(htep) + Math.Sqrt(hrep))), 2.0); //Equation A.2.4
                    adft = aem;
                    //Use A.3 to calculate Ldft
                    Ldft = a3_loss(freq, adft, Tpol, htep, hrep, omega, dn);  //Subrountine for A.3
                    if (Ldft < 0)
                        Ldsph = 0;
                    else
                        Ldsph = (1 - (hsph / hrep)) * Ldft; //Equation A.2.5
                }
            }
            return Ldsph;
        }

        //P.2001-2-2015: Attachement A.3: First-Term Spherical Earth Diffraction Loss using A.3 to find Ldft
        double a3_loss(double freq, double adft, double Tpol, double htep, double hrep, double omega, double dn)
        {
            const double erland = 22.0;  //Relative Permittivity for Land
            const double ersea = 80.0; //Relative Permittivity fo Sea
            const double sigmaland = 0.003; //Conducitivty for land (S/m)
            const double sigmasea = 5.0; //Conductivity for sea (S/m)
            double Ldft, Ldftland, Ldftsea;

            // LAND
            Ldftland = a3_subLdft(adft, freq, erland, sigmaland, Tpol, htep, hrep, dn); //Calculate Ldft using A.3.2 and A.3.8 to get Ldftland
            // SEA
            Ldftsea = a3_subLdft(adft, freq, ersea, sigmasea, Tpol, htep, hrep, dn); //Calculate Ldft using A.3.2 and A.3.8 to get Ldftsea
            Ldft = omega * Ldftsea + (1 - omega) * Ldftland; //Equation A.3.1
            return Ldft;
        }

        //P.2001-2-2015: Attachement A.3:Equations A.3.2-A.3.8 First-Term Spherical Earth Diffraction Loss using A.3 to find Ldft
        double a3_subLdft(double adft, double freq, double er, double sigma, double Tpol, double htep, double hrep, double dn)
        {
            double subLdft, KH, KV, Beta, X, Yt, Yr, FX;
            double Bt, Br, GYt, GYr;
            KH = 0.036 * Math.Pow(adft * freq, -1.0 / 3.0) * Math.Pow(Math.Pow(er - 1, 2.0) + Math.Pow(18 * sigma / freq, 2.0), -1.0 / 4.0); //Equation A.3.2a
            KV = KH * Math.Pow(Math.Pow(er, 2.0) + Math.Pow(18 * sigma / freq, 2.0), 1.0 / 2.0); //Equation A.3.2b
            if (Tpol == 1)
                Beta = (1.0 + 1.6 * Math.Pow(KV, 2.0) + 0.67 * Math.Pow(KV, 4.0)) / (1.0 + 4.5 * Math.Pow(KV, 2.0) + 1.53 * Math.Pow(KV, 4.0)); //Equation A.3.3
            else
                Beta = (1 + 1.6 * Math.Pow(KH, 2.0) + 0.67 * Math.Pow(KH, 4.0)) / (1 + 4.5 * Math.Pow(KH, 2.0) + 1.53 * Math.Pow(KH, 4.0)); //Equation A.3.3
            X = 21.88 * Beta * Math.Pow((freq) / Math.Pow(adft, 2.0), 1.0 / 3.0) * dn; //Equation A.3.4
            Yt = 0.9575 * Beta * Math.Pow(Math.Pow(freq, 2.0) / adft, 1.0 / 3.0) * htep; //Equation A.3.5a
            Yr = 0.9575 * Beta * Math.Pow(Math.Pow(freq, 2.0) / adft, 1.0 / 3.0) * hrep; //Equation A.3.5b
            if (X >= 1.6)
                FX = 11 + 10 * Math.Log10(X) - 17.6 * X; //Equation A.3.6
            else
                FX = -20 * Math.Log10(X) - 5.6488 * Math.Pow(X, 1.425); //Equation A.3.6
            Bt = Beta * Yt; //Equation A.3.7a
            if (Bt > 2.0)
                GYt = 17.6 * Math.Pow(Bt - 1.1, 0.5) - 5 * Math.Log10(Bt - 1.1) - 8; //Equation A.3.7
            else
                GYt = 20 * Math.Log10(Bt + 0.1 * Math.Pow(Bt, 3.0)); //Equation A.3.7
            Br = Beta * Yr; //Equation A.3.7a
            if (Br > 2)
                GYr = 17.6 * Math.Pow(Br - 1.1, 0.5) - 5 * Math.Log10(Br - 1.1) - 8; //Equation A.3.7
            else
                GYr = 20 * Math.Log10(Br + 0.1 * Math.Pow(Br, 3.0)); //Equation A.3.7
            //Limit G(Y) such that G(Y)>= 2+20*log10(K)
            if (Tpol == 1) //Vertical
            {
                if (GYt < (2 + 20 * Math.Log10(KV)))
                    GYt = 2 + 20 * Math.Log10(KV);
                if (GYr < (2 + 20 * Math.Log10(KV)))
                    GYr = 2 + 20 * Math.Log10(KV);
            }
            else //Horizontal
            {
                if (GYt < (2 + 20 * Math.Log10(KH)))
                    GYt = 2 + 20 * Math.Log10(KH);
                if (GYr < (2 + 20 * Math.Log10(KH)))
                    GYr = 2 + 20 * Math.Log10(KH);
            }
            subLdft = -FX - GYt - GYr; //Equation A.3.8 %Output
            return subLdft;
        }

        // P.2001-2-2015: Attachement B.5: Percentage of time a given clear-air fade level is exceeded on a troposcatter path
        double b5_Qcaftropo(double A)
        {
            double Qcaftropo;
            if (A < 0.0)
                Qcaftropo = 100.0; // Equation B.5.1a
            else
                Qcaftropo = 0.0; // Equation B.5.1b
            return Qcaftropo;
        }

        double Qiter_sm3(double A, double Q0ra, double Pr6, double[] Pm, double[] Gm, double hRtop, double hrainlo,
            double dr, double kmod, double alphamod, double a1, double b1, double c1)
        {
            double Qcaftropo, Qrain, ret;
            Qcaftropo = b5_Qcaftropo(A);
            Qrain = c3_Qrain(A, hrainlo, hRtop, Pm, Gm, Pr6, dr, kmod, alphamod, a1, b1, c1); // C.3: Calculate Qrain(A)
            ret = Qrain * (Q0ra / 100.0) + Qcaftropo * (1.0 - (Q0ra / 100.0));  // Qiter(A):Equation 4.1.3    
            return ret;
        }

        // P.2001-2-2015: Attachement B.3: Calculation of the Notional Zero-Fade Annual Percentage Time
        double b3_Q0ca(double dca, double eca, double hca, double freq, double phimn, double K)
        {
            double Q0ca, qw, Cg;
            qw = K * (Math.Pow(dca, 3.1)) * (Math.Pow(1 + eca, -1.29)) * (Math.Pow(freq, 0.8)) * (Math.Pow(10.0, -0.00089 * hca)); // Equation B.3.1
            if (Math.Abs(phimn) <= 45.0)
                Cg = 10.5 - (5.6 * Math.Log10(1.1 + Math.Pow(Math.Abs(Math.Cos(deg2rad(2 * phimn))), 0.7))) -
                    (2.7 * Math.Log10(dca)) + (1.7 * Math.Log10(1.0 + eca)); // Equation B.3.2a
            else
                Cg = 10.5 - (5.6 * Math.Log10(1.1 - Math.Pow(Math.Abs(Math.Cos(deg2rad(2 * phimn))), 0.7))) -
                    (2.7 * Math.Log10(dca)) + (1.7 * Math.Log10(1.0 + eca)); // Equation B.3.2b
            if (Cg > 10.8)
                Cg = 10.8;
            Q0ca = Math.Pow(10.0, -0.1 * Cg) * qw; // Equation B.3.3
            return Q0ca;
        }

        double deg2rad(double deg)
        {
            const double d2r = Math.PI / 180.0;
            return deg * d2r;
        }

        public double get_txt_value(double lat, double lon, double[,] temp)
        {
            int rptr, cptr;
            double cv, rv;
            int rdim = temp.GetUpperBound(0), cdim = temp.GetUpperBound(1);
            if (cdim == 719) // For Troposcatter, no bi-linear interpolation
            {
                cv = (double)cdim * (lon + 179.75) / 359.5;
                rv = (double)rdim * (lat + 89.75) / 179.5;
                rptr = (int)(rdim - Math.Round(rv));
                cptr = (int)Math.Round(cv);
                return temp[rptr, cptr];
            }
            else // interpolate
            {
                if (lon < 0.0)
                    cv = (lon + 360.0) / (360.0 / (double)cdim);
                else
                    cv = lon / (360.0 / (double)cdim);
                rv = (double)rdim - (90.0 + lat) / (180.0 / (double)rdim);
                int r1 = (int)Math.Floor(rv), r2 = (int)Math.Ceiling(rv);
                int c1 = (int)Math.Floor(cv), c2 = (int)Math.Ceiling(cv);
                if (r1 == r2 && c1 == c2)
                    return temp[r1, c1];
                else if (r1 == r2 && c1 != c2)
                    return (temp[r1, c2] - temp[r1, c1]) * (cv - c1) + temp[r1, c1];
                else if (r1 != r2 && c1 == c2)
                    return (temp[r2, c1] - temp[r1, c1]) * (rv - r1) + temp[r1, c1];
                else
                {
                    double sc1 = (temp[r1, c2] - temp[r1, c1]) * (cv - c1) + temp[r1, c1];
                    double sc2 = (temp[r2, c2] - temp[r2, c1]) * (cv - c1) + temp[r2, c1];
                    double sr = (sc2 - sc1) * (rv - r1) + sc1;
                    return sr;
                }
            }
        }

        void c2_calc_mergeC5(double phin, double phie, double q_time, double hrainlo, double hrainhi, double drain, double Tpol, double freq, double hlo, double hhi,
            out double Fwvr, out double Q0ra, out double[] Pm, out double[] Gm, out double hRtop, out double Pr6, out double dr, out double kmod,
            out double alphamod, out double a1, out double b1, out double c1)
        {
            //Get Pr6 from 'Esarain_Pr6_v5.txt' using phin and phie
            //Get MT from 'Esarain_MT_v5.txt' using phin and phie
            //Get Brain from 'Esarain_Beta_v5.txt' using phin and phie
            double Qtran, Q = 0.0;

            Pr6 = get_txt_value(phin, phie, Esarain_Pr6_v5);
            double MT = get_txt_value(phin, phie, Esarain_Mt_v5);
            double Brain = get_txt_value(phin, phie, Esarain_Beta_v5);
            double h0 = get_txt_value(phin, phie, h0text);
            double hR = 360 + 1000 * h0; //Equation C.2.1
            hRtop = hR + 2400; //Equation C.2.2

            if (Pr6 == 0 || hrainlo >= hRtop) //No Rain,
            {
                Q0ra = 0; //Main outputs
                Fwvr = 0; //Main outputs
                Pm = new double[] { 0 }; //Place Holders when there is no rain
                Gm = new double[] { 0 }; //Place Holders when there is no rain
                dr = 0; //Place Holders when there is no rain
                kmod = 0; //Place Holders when there is no rain
                alphamod = 0; //Place Holders when there is no rain
                a1 = 0; //Place Holders when there is no rain
                b1 = 0; //Place Holders when there is no rain
                c1 = 0; //Place Holders when there is no rain
            }
            else
            {
                double Mc, MS, tau, erain;
                double alpha, k, k1ghz, alpha1ghz, drmin;
                double kmod1, kmod2, alphamod1, alphamod2, alphamod3;
                Mc = Brain * MT; //Equation C.2.3a {Correct}
                MS = (1.0 - Brain) * MT; //Equation C.2.3b
                Q0ra = Pr6 * (1.0 - Math.Exp((-0.0079 * MS) / Pr6)); //Equation C.2.4
                a1 = 1.09; //C.2.5a
                b1 = (Mc + MS) / (21797 * Q0ra); //Equation C.2.5b
                c1 = 26.02 * b1; //Equation C.2.5c
                Qtran = Q0ra * Math.Exp((a1 * (2 * b1 - c1)) / (c1 * c1)); //Equation C.2.6
                if (Tpol == 1) //Vertical
                    tau = 90.0;
                else
                    tau = 0.0;
                erain = (0.001 * (hrainhi - hrainlo)) / drain; //Equation C.2.7
                //Use ITU-R P.838 to calculate rain regression coefficients for k and alpha
                if (freq < 1.0)
                {
                    rain_coefficients_838(1, deg2rad(tau), erain, out k1ghz, out alpha1ghz);
                    k = freq * k1ghz; //Equation C.2.8a //k1ghz are calculated for 1GHz
                    alpha = alpha1ghz; //Equation C.2.8b //alpha1ghz are calculated for 1GHz
                }
                else
                    rain_coefficients_838(freq, deg2rad(tau), erain, out k, out alpha);

                dr = Math.Min(drain, 300.0); //Equation C.2.9a
                drmin = Math.Max(dr, 1.0); //Equation C.2.9b
                kmod1 = Math.Pow(1.763, (alpha)) * k;
                kmod2 = (0.6546 * Math.Exp(-0.009516 * drmin) + 0.3499 * Math.Exp(-0.001182 * drmin));
                kmod = kmod1 * kmod2; //Equation C.2.10a
                alphamod1 = (0.753 + (0.197 / drmin)) * alpha;
                alphamod2 = 0.1572 * Math.Exp(-0.02268 * drmin);
                alphamod3 = 0.1594 * Math.Exp(-0.0003617 * drmin);
                alphamod = alphamod1 + alphamod2 - alphamod3; //Equation C.2.10b
                Pm = new double[49]; //zeros(1,49); //Proability of Particular case, initialize all Pm to zero
                Gm = new double[49]; //zeros(1,49);
                Gm[0] = 1.0;
                double[] Hn = new double[49];
                int m, dvptr = 0;
                for (double dv = -2400.0; dv <= 2400.0; dv += 100.0)
                    Hn[dvptr++] = dv;
                m = 0;  //marker
                //For each line of Table C.2.1, for n from 1 to 49, do the following:
                for (int n = 0; n < Pm.Length; ++n)
                {
                    //a)
                    //clear hT; //Equation C2. a)
                    double hT, deltah = 0.0, GAMMAslice, G, Pm_temp;
                    int slo, shi;
                    hT = hR + Hn[n]; //Hn is the corresponding relative height entry in Table C.2.1
                    //b)
                    if (hrainlo >= hT) //Equation C2. b)
                    {
                        //repeat from a) for the next value of n
                        //Do nothing
                    }
                    else
                    {
                        //do c)
                        if (hrainhi > (hT - 1200)) //Equation C2. c)
                        {
                            //c.i)Use the method C.5 to set Gm to the path-averaged multiplier for this path geometry relative ot the melting layer
                            //clear slo;
                            //clear shi;
                            slo = 1 + Convert.ToInt32(Math.Floor((hT - hlo) / 100)); //Equation C.5.1a
                            shi = 1 + Convert.ToInt32(Math.Floor((hT - hhi) / 100)); //Equation C.5.1b
                            if (slo < 1.0)
                                G = 0;
                            //End of G calculation
                            else if (shi > 12.0)
                                G = 1;
                            //End of G claculation
                            else if (slo == shi)
                            {
                                //clear deltah;
                                deltah = 0.5 * (hlo + hhi) - hT; //Equation C.5.2
                                GAMMAslice = c4_Gamma(deltah);  //Equation C.5.2
                                G = GAMMAslice;
                            }   //End of claculation
                            else
                            {
                                int[] s = new int[2];
                                G = 0.0; //initialize, Equation C.5.3
                                s[0] = Math.Max(shi, 1); //Equation C.5.4a //sfirst
                                s[1] = Math.Min(slo, 12); //Equation C.5.4b //slast
                                for (int j = s[0]; j <= s[1]; ++j)
                                {
                                    //Condition 1
                                    if (shi < j && j < slo)
                                    {
                                        deltah = 100.0 * (0.5 - j); //Equation C.5.5a
                                        Q = 100.0 / (hhi - hlo); //Equation C.5.5b
                                    }
                                    else if (j == slo) //Condition 2
                                    {
                                        deltah = 0.5 * ((hlo - hT - (100 * (j - 1)))); //Equation C.5.6a
                                        Q = (hT - (100 * (j - 1)) - hlo) / (hhi - hlo); //Equation C.5.6b
                                    }
                                    else if (j == shi) //Condition 3
                                    {
                                        deltah = 0.5 * (hhi - hT - (100 * j)); //Equation C.5.7a
                                        Q = (hhi - (hT - 100 * j)) / (hhi - hlo); //Equation C.5.7b
                                    }
                                    GAMMAslice = c4_Gamma(deltah); //Equation C.5.8;
                                    G = G + (Q * GAMMAslice); //Equation C.5.9
                                }
                            }
                            if (slo > 12)
                            {
                                Q = (hT - 1200 - hlo) / (hhi - hlo); //Equation C.5.10
                                G = G + Q; //Equation C.5.11
                            }
                            Gm[m] = G;
                            //c.ii) set Pm= to the third column (Cap pi sub n) in table C.2.1
                            Pm[m] = prob_C21(n);
                            //c.iii)
                            if (n < 49)
                                m = m + 1; //add 1 to array index m
                            //c.iv) repeat from a) for the next value of n
                        }
                        else
                        {
                            //Continue to d)
                            //d) Accumulate (Cap pi sub n) from Table C.2.1 into Pm,
                            //   set Gm=1 and repeat a) for the next value of n.
                            if (Pm[m] == 0.0)
                                Pm[m] = prob_C21(n);
                            else
                            {
                                Pm_temp = Pm[m];
                                Pm[m] = Pm_temp + prob_C21(n);
                            }
                            Gm[m] = 1.0;
                        }
                    }
                }

                //Set the number of values in arrays Gm and Pm according to M=m Equation C.2.12
                Array.Resize(ref Pm, m + 1);  //clear temp_Pm; temp_Pm=Pm(1:m); clear Pm; Pm=temp_Pm;
                Array.Resize(ref Gm, m + 1);  //clear temp_Gm; temp_Gm=Gm(1:m); clear Gm; Gm=temp_Gm;
                double Rwvr = 6.0 * ((Math.Log10(Q0ra / q_time)) / (Math.Log10(Q0ra / Qtran))) - 3; //Equation C.2.13a
                double temp_sum = 0.0;
                for (int i = 0; i < Gm.Length; ++i)
                    temp_sum += Gm[i] * Pm[i];//Equation C.2.13pre
                Fwvr = 0.5 * (1.0 + Math.Tanh(Rwvr)) * temp_sum; //Equation C.2.13 {Correct}
            }
        }

        /// <summary>
        /// ITU-R P.838-3 Rain Attenuation
        /// </summary>
        void rain_coefficients_838(double freq, double tau, double erain, out double k, out double alpha)
        {
            // Use ITU-R P.838 to calculate rain regression coefficients for k and alpha
            // Tables from ITU-R 838.3 (Tables 1-4)
            double[] kH_aj = { -5.33980, -0.35351, -0.23789, -0.94158 };
            double[] kH_bj = { -0.10008, 1.26970, 0.86036, 0.64552 };
            double[] kH_cj = { 1.13098, 0.45400, 0.15354, 0.16817 };
            double kH_mk = -0.18961;
            double kH_ck = 0.71147;

            double[] kV_aj = { -3.80595, -3.44965, -0.39902, 0.50167 };
            double[] kV_bj = { 0.56934, -0.22911, 0.73042, 1.07319 };
            double[] kV_cj = { 0.81061, 0.51059, 0.11899, 0.27195 };
            double kV_mk = -0.16398;
            double kV_ck = 0.63297;

            double[] aH_aj = { -0.14318, 0.29591, 0.32177, -5.37610, 16.1721 };
            double[] aH_bj = { 1.82442, 0.77564, 0.63773, -0.96230, -3.29980 };
            double[] aH_cj = { -0.55187, 0.19822, 0.13164, 1.47828, 3.43990 };
            double aH_ma = 0.67849;
            double aH_ca = -1.95537;

            double[] aV_aj = { -0.07771, 0.56727, -0.20238, -48.2991, 48.5833 };
            double[] aV_bj = { 2.33840, 0.95545, 1.14520, 0.791669, 0.791459 };
            double[] aV_cj = { -0.76284, 0.54039, 0.26809, 0.116226, 0.116479 };
            double aV_ma = -0.053739;
            double aV_ca = 0.83433;

            // clear temp_kH;
            double temp_kH = 0.0;
            for (int j = 0; j < 4; ++j)
                temp_kH += kH_aj[j] * Math.Exp(-1.0 * Math.Pow((Math.Log10(freq) - kH_bj[j]) / kH_cj[j], 2.0));  //ITU-R 838.3 Equation 2
            double sum_kH = temp_kH + kH_mk * Math.Log10(freq) + kH_ck; //ITU-R 838.3 Equation 2
            double kH = Math.Pow(10, sum_kH); //ITU-R 838.3 Equation 2

            double temp_alphaH = 0.0;
            for (int j = 0; j < 5; ++j)
                temp_alphaH += aH_aj[j] * Math.Exp(-1.0 * Math.Pow((Math.Log10(freq) - aH_bj[j]) / aH_cj[j], 2.0));  //ITU-R 838.3 Equation 3
            double alphaH = temp_alphaH + aH_ma * Math.Log10(freq) + aH_ca; //ITU-R 838.3 Equation 3

            double temp_kV = 0.0;
            for (int j = 0; j < 4; ++j)
                temp_kV += kV_aj[j] * Math.Exp(-1.0 * Math.Pow((Math.Log10(freq) - kV_bj[j]) / kV_cj[j], 2.0));  //ITU-R 838.3 Equation 2
            double sum_kV = temp_kV + kV_mk * Math.Log10(freq) + kV_ck; //ITU-R 838.3 Equation 2
            double kV = Math.Pow(10.0, sum_kV); // ITU-R 838.3 Equation 2

            double temp_alphaV = 0.0;
            for (int j = 0; j < 5; ++j)
                temp_alphaV += aV_aj[j] * Math.Exp(-1.0 * Math.Pow((Math.Log10(freq) - aV_bj[j]) / aV_cj[j], 2.0));  //ITU-R 838.3 Equation 3
            double alphaV = temp_alphaV + aV_ma * Math.Log10(freq) + aV_ca; //ITU-R 838.3 Equation 3

            k = (kH + kV + (kH - kV) * Math.Pow(Math.Cos(erain), 2) * Math.Cos(2 * tau)) / 2; //ITU-R 838.3 Equation 4
            alpha = (kH * alphaH + kV * alphaV + (kH * alphaH - kV * alphaV) * Math.Pow(Math.Cos(erain), 2.0) * Math.Cos(2 * tau)) / (2 * k); //ITU-R 838.3 Equation 5
        }

        /// <summary>
        /// P.2001-2-2015: Find Value of Table C.2.1
        /// </summary>
        /// <param name="index">Zero Based Index</param>
        /// <returns>Probability Based on Table C.2.1</returns>
        double prob_C21(int index)
        {
            return probability_C21[index];
        }


        /// <summary>
        /// P.2001-2-2015: Attachement C.3: Percentage Time a Given Precipitation Fade Level is Exceeded
        /// </summary>
        double c3_Qrain(double A, double hrainlo, double hRtop, double[] Pm, double[] Gm, double Pr6,
            double dr, double kmod, double alphamod, double a1, double b1, double c1)
        {
            double Qrain = 0.0;
            if (A < 0)
                Qrain = 100.0; //C.3.1a  
            else if (A >= 0) //This is dependent on whether the path is classifed as "rain" or "non rain"  
            {
                if (Pr6 == 0 || hrainlo >= hRtop) //No Rain, %Pulled from C.2
                    Qrain = 0; //Equation C.3.1b, "Non-rain"
                else
                {
                    //"Rain", Where these Pm and Gm are defined in C.2.11-C.2.13a
                    double drlim = Math.Max(dr, 0.001);  //Equation C.3.1e [km]
                    double temp_Qrain = 0.0;
                    double Rm;
                    for (int m = 0; m < Gm.Length; ++m)
                    {
                        Rm = Math.Pow(A / (Gm[m] * drlim * kmod), 1.0 / alphamod); //Equation C.3.1d
                        temp_Qrain += (Pm[m] * Math.Exp(-1 * (a1 * Rm * (b1 * Rm + 1.0)) / (c1 * Rm + 1.0))); //Equation C.3.1c
                    }
                    Qrain = 100.0 * temp_Qrain;
                }
            }
            return Qrain;
        }

        /// <summary>
        /// P.2001-2-2015: C.4: Melting Layer Model, returning the attenuation multiplier
        /// </summary>
        double c4_Gamma(double deltah)
        {
            double capGamma = 0.0;
            if (0.0 < deltah)
                capGamma = 0; //Equation C.4.1
            else if (-1200.0 <= deltah && deltah <= 0.0)
            {
                double gamma1 = 4 * Math.Pow((1.0 - Math.Exp(deltah / 70)), 2.0);
                double gamma2 = Math.Pow((1.0 - Math.Exp(-(Math.Pow((deltah / 600), 2.0)))), 2.0);
                double gamma3 = (4 * (Math.Pow((1.0 - Math.Exp(deltah / 70)), 2.0)) - 1);
                capGamma = gamma1 / (1 + (gamma2 * gamma3));
            }
            else if (deltah < -1200.0)
                capGamma = 1; //Equation C.4.1
            return capGamma;
        }

        // P.2001: D.1: Find the Distance for dlm/dtm
        void find_ITU_dist(double rxlat, double rxlon, double txlat, double txlon, double[] h, double dn,
            out double dtm, out double dlm, out double dct, out double dcr, out double omega)
        {
            int mid_npts_rx, mid_npts_tx;
            int marker1, ind1, ind2;
            double consec_len_rx, consec_len_tx;
            double[] tropo_value_rx = new double[h.Length], tropo_value_tx = new double[h.Length];
            double[] r_lat = new double[h.Length], r_lon = new double[h.Length], d_tmp = new double[h.Length];
            double[] t_lat = new double[h.Length], t_lon = new double[h.Length], t_tmp = new double[h.Length];

            List<double> consec_dist_rx = new List<double>(), conesc_dist_tx = new List<double>();
            Vincenty vin = new Vincenty();
            GenerateCoords(rxlat, rxlon, txlat, txlon, h.Length, out r_lat, out r_lon, out d_tmp);
            mid_npts_rx = h.Length;// This is based upon the number of points in the terrain model (h).
            for (int i = 0; i < r_lat.Length; ++i)
                tropo_value_rx[i] = get_txt_value(r_lat[i], r_lon[i], TropoClim); //Gets the Climate of each point starting at the rx
            marker1 = 0;
            ind1 = 0;
            ind2 = 0;
            for (int i = 0; i < tropo_value_rx.Length; ++i)
            {
                if (tropo_value_rx[i] != 0.0)
                    ind2 = i;
                else
                {
                    if (ind2 == 0)
                        ind2 = i;
                    else
                    {
                        consec_dist_rx.Add(vin.GetNewDist_m_SphericalEarth(r_lat[ind1], r_lon[ind1], r_lat[ind2], r_lon[ind2]) / 1000.0);
                        ++marker1;
                        ind1 = i;
                        if (i < tropo_value_rx.Length - 1)
                            ind1 = i + 1;
                    }
                }
            }
            consec_dist_rx.Add(vin.GetNewDist_m_SphericalEarth(r_lat[ind1], r_lon[ind1], r_lat[ind2], r_lon[ind2]) / 1000.0);
            consec_len_rx = consec_dist_rx.Max();
            dcr = consec_dist_rx[0]; //The r_lat/lon starts at the receiver
            mid_npts_tx = h.Length;// This is based upon the number of points in the terrain model (h).
            GenerateCoords(txlat, txlon, rxlat, rxlon, h.Length, out t_lat, out t_lon, out d_tmp);
            for (int i = 0; i < t_lat.Length; ++i)
                tropo_value_tx[i] = get_txt_value(t_lat[i], t_lon[i], TropoClim); //Gets the Climate of each point starting at the rx
            marker1 = 0;
            ind1 = 0;
            ind2 = 0;
            for (int i = 0; i < tropo_value_tx.Length; ++i)
            {
                if (tropo_value_tx[i] != 0.0)
                    ind2 = i;
                else
                {
                    if (ind2 == 0)
                        ind2 = i;
                    else
                    {
                        conesc_dist_tx.Add(vin.GetNewDist_m_SphericalEarth(t_lat[ind1], t_lon[ind1], t_lat[ind2], t_lon[ind2]) / 1000.0);
                        ++marker1;
                        ind1 = i;
                        if (i < tropo_value_tx.Length - 1)
                            ind1 = i + 1;
                    }
                }
            }

            conesc_dist_tx.Add(vin.GetNewDist_m_SphericalEarth(t_lat[ind1], t_lon[ind1], t_lat[ind2], t_lon[ind2]) / 1000.0);
            consec_len_tx = conesc_dist_tx.Max();
            dct = conesc_dist_tx[0]; //The t_lat/lon starts at the receiver

            //Need to refine how to define dtm, dlm based on h.

            //Find the first inland point near the sea, and go inland, checking the
            //elevation to see when it become greater than 100m, and then see how long that path if,
            //checking that it does not exceed 50km.
            dtm = consec_len_tx;  //Coast and inland distance.
            dlm = consec_len_tx;  //Inland A2, greater than 100m, not to exceed 50km inland.
            omega = Math.Round(1.0 - (consec_len_tx / dn), 2);
        }

    }
}
