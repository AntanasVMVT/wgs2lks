using System;


//Original source - https://github.com/vilnius/lks2wgs/
public class GeoConverter {
	private const double Pi = Math.PI;
	private const double A = 6378137;
	private const double B = A * (1 - (1 / 298.257223563));
	private const double C = (A * A - B * B) / (A * A);
	private const double D = (A - B) / (A + B);
	private const double E = (1 - (C / 4) - (3 * C * C / 64) - (5 * C * C * C / 256));
	private const double F = ((3.0 / 8) * (C + (C * C / 4) + (15 * C * C * C / 128)));
	private const double G = ((15.0 / 256) * (C * C + (3 * C * C * C / 4)));
	private const double H = (180 * (A * (1 - D) * (1 - D * D) * (1 + (9.0 / 4) * D * D + (255.0 / 64) * D * D * D * D) * (Pi / 180)));
	private const double I = (35 * C * C * C / 3072);
	private const double J = A * (1 - C);
	private const double K = 0.9998;
	private const double L = ((3 * D / 2) - (27 * D * D * D / 32));
	private const double M = ((21 * D * D / 16) - (55 * D * D * D * D / 32));
	private const double N = (151 * D * D * D / 96);
	private const double O = (1097 * D * D * D * D / 512);

	// WGS-84 to LKS-94
	public static (double east, double north) GeoToGrid(double lat, double lon, int round = 2) {
		double latrad = DegToRad(lat);
		double lonrad = DegToRad(lon - 24);
		double coslat = Math.Cos(latrad);
		double sinlat = Math.Sin(latrad);
		double lr2c = 1 - (C * sinlat * sinlat);
		double t = Math.Tan(latrad); var t2 = t * t; var t4 = Math.Pow(t, 4); var t6 = Math.Pow(t, 6);
		double nu = A / Math.Sqrt(lr2c); var nusin = nu * sinlat;
		double psi = nu / J / Math.Pow(lr2c, 1.5); var psi2 = psi * psi;
		double m = A * (E * latrad - F * Math.Sin(2 * latrad) + G * Math.Sin(4 * latrad) - I * Math.Sin(6 * latrad));

		return (
			Math.Round(500000 + (K * nu * lonrad * coslat * (1 + ((Math.Pow(lonrad, 2) / 6) * coslat * coslat * (psi - t2)) +
				((Math.Pow(lonrad, 4) / 120) * Math.Pow(coslat, 4) * (4 * Math.Pow(psi, 3) * (1 - 6 * t2) + psi2 * (1 + 8 * t2) - psi * 2 * t2 + t4)) +
				((Math.Pow(lonrad, 6) / 5040) * Math.Pow(coslat, 6) * (61 - 479 * t2 + 179 * t4 - t6)))), round),
			Math.Round(K * (m + ((Math.Pow(lonrad, 2) / 2) * nusin * coslat) + ((Math.Pow(lonrad, 4) / 24) * nusin * Math.Pow(coslat, 3) * (4 * psi2 + psi - t2)) +
				((Math.Pow(lonrad, 6) / 720) * nusin * Math.Pow(coslat, 5) * (8 * Math.Pow(psi, 4) * (11 - 24 * t2) - 28 * Math.Pow(psi, 3) * (1 - 6 * t2) + psi2 * (1 - 32 * t2) - psi * 2 * t2 + t4)) +
				((Math.Pow(lonrad, 8) / 40320) * nusin * Math.Pow(coslat, 7) * (1385 - 3111 * t2 + 543 * t4 - t6))), round)
			);
	}

	// LKS-94 to WGS-84
	public static (double lat, double lon) GridToGeo(double x_east, double y_north, int round = 6) {
		double east = x_east - 500000;
		double sigma = ((y_north / K) * Pi) / H;
		double footlat = sigma + L * Math.Sin(2 * sigma) + M * Math.Sin(4 * sigma) + N * Math.Sin(6 * sigma) + O * Math.Sin(8 * sigma);
		double fl2c = 1 - (C * Math.Pow(Math.Sin(footlat), 2));
		double rho = J / Math.Pow(fl2c, 1.5);
		double nu = A / Math.Sqrt(fl2c);
		double psi = nu / rho; var psi2 = psi * psi;
		double t = Math.Tan(footlat); var t2 = t * t; var t4 = Math.Pow(t, 4); var t6 = Math.Pow(t, 6);
		double xval = east / (K * nu);
		double rhokt = t / (K * rho);
		double secLat = 1 / Math.Cos(footlat);

		return (
			Math.Round(RadToDeg(footlat - (rhokt * (east * xval / 2)) + (rhokt * (east * Math.Pow(xval, 3) / 24) * (-4 * psi2 + 9 * psi * (1 - t2) + 12 * t2)) -
				(rhokt * (east * Math.Pow(xval, 5) / 720) * (8 * Math.Pow(psi, 4) * (11 - 24 * t2) - 12 * Math.Pow(psi, 3) * (21 - 71 * t2) + 15 * psi2 * (15 - 98 * t2 + 15 * t4) + 180 * psi * (5 * t2 - 3 * t4) + 360 * t4)) +
				(rhokt * (east * Math.Pow(xval, 7) / 40320) * (1385 + 3633 * t2 + 4095 * t4 + 1575 * t6))), round),
			Math.Round(RadToDeg(DegToRad(24) + ((xval * secLat) - ((Math.Pow(xval, 3) / 6) * secLat * (psi + 2 * t2)) +
				((Math.Pow(xval, 5) / 120) * secLat * (-4 * Math.Pow(psi, 3) * (1 - 6 * t2) + psi2 * (9 - 68 * t2) + 72 * psi * t2 + 24 * t4)) -
				((Math.Pow(xval, 7) / 5040) * secLat * (61 + 662 * t2 + 1320 * t4 + 720 * t6)))), round)
			);
	}

	public static double DegToRad(double deg) => deg * (Pi / 180);
	public static double RadToDeg(double rad) => rad * (180 / Pi);
}
