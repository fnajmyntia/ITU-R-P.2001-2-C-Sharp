# ITU-R-P.2001-2-C-Sharp

This is a C# implementation of the ITU-R P.2001-2,2015 RF Propagation Model

The main routine is:
  public double p2001(double freq, double Tpol, double txlat, double txlon, double rxlat, double rxlon,
         double Hrg, double Htg, double Tpc, double Gr, double Gt, double[] h)
Parameters:
       freq: (GHz) 30 MHz - 50 GHz
       Tpol: Antenna Polarization: Vertical (1), Horizontal (0)
       txlat: (degrees) Latitude of Transmitter
       txlon: (degrees) Longitude of Transmitter
       rxlat: (degrees) Latitude of Receiver
       rxlon: (degrees) Longitude of Receiver
       Hrg: (meters) Height of Electrical Center of Receiving Antenna above ground
       Htg: (meters) Height of Electrical Center of Transmitting Antenna above ground
       Tpc: Percentage of average year for which the prdicted basic transmission loss is not exceeded (%), Example 0.00001% to 99.99999%
       Gr: Gain (dBi) of receiving antenna in the azimuthal direction of the path towards the other antenna, and at the elevation angle about the local horizontal of the other antenna in the case of a line-of-sight path, otherwise of the antenna's radio horizon, for median effective Earth radius.
       Gt: Gain (dBi) of transmitting antenna in the azimuthal direction of the path towards the other antenna, and at the elevation angle about the local horizontal of the other antenna in the case of a line-of-sight path, otherwise of the antenna's radio horizon, for median effective Earth radius.
       h: (meters) Elevation points, evenly spaced.

Requires data files:
DN_Median.csv
DN_SubSlope.csv
DN_SupSlope.csv
Esarain_Beta_v5.csv
Esarain_Mt_v5.csv
Esarain_Pr6_v5.csv
FoEs0.1.csv
FoEs01.csv
FoEs10.csv
FoEs50.csv
TropoClim.csv
dndz_01.csv
h0.csv
surfwv_50_fixed.csv

For ITU sotware and validation examples see https://www.itu.int/en/ITU-R/study-groups/rsg3/Pages/iono-tropo-spheric.aspx
For publication: Recommendation P.2001-2 (07/2015) see https://www.itu.int/rec/R-REC-P.2001-2-201507-I/en
