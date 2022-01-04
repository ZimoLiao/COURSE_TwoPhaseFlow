#include "LevelSetSolver.h"

void LevelSetSolver::QUICK(Array& phi)
{
	double aw, ae, as, an;
	double phiP, phiWW, phiW, phiE, phiEE, phiSS, phiS, phiN, phiNN;
	double fij;
	for (int i = 0; i != n_; i++) {
		for (int j = 0; j != n_; j++) {
			aw = double(u.w(i, j) > 0);
			ae = double(u.e(i, j) > 0);
			as = double(u.s(i, j) > 0);
			an = double(u.n(i, j) > 0);

			fij = 0.0;
			phiP = phi(i, j);
			phiWW = phi(i - 2, j);
			phiW = phi(i - 1, j);
			phiE = phi(i + 1, j);
			phiEE = phi(i + 2, j);
			phiSS = phi(i, j - 2);
			phiS = phi(i, j - 1);
			phiN = phi(i, j + 1);
			phiNN = phi(i, j + 2);

			fij += 0.125 * u.e(i, j) * (-ae * phiW + (3.0 + 3.0 * ae) * phiP\
				+ (6.0 - 3.0 * ae) * phiE - (1.0 - ae) * phiEE);
			fij -= 0.125 * u.w(i, j) * (-aw * phiWW + (3.0 + 3.0 * aw) * phiW\
				+ (6.0 - 3.0 * aw) * phiP - (1.0 - aw) * phiE);
			fij += 0.125 * u.n(i, j) * (-an * phiS + (3.0 + 3.0 * an) * phiP\
				+ (6.0 - 3.0 * an) * phiN - (1.0 - an) * phiNN);
			fij -= 0.125 * u.s(i, j) * (-as * phiSS + (3.0 + 3.0 * as) * phiS\
				+ (6.0 - 3.0 * as) * phiP - (1.0 - as) * phiN);

			fij *= -dt_ / h_;

			f(i, j) = fij;
		}
	}
}

void LevelSetSolver::HJENO(Array& phi)
{
}

void LevelSetSolver::F(Array& phi)
{
	switch (scheme_s_)
	{
	case 1:
		QUICK(phi);
		break;
	case 2:
		HJENO(phi);
		break;
	default:
		QUICK(phi);
		break;
	}
}

void LevelSetSolver::EULER()
{
	F(phin);
	phin += f;

	phin.ExtrapGhost(0);
}

void LevelSetSolver::TVDRK2()
{
	phin1 = phin;
	F(phin);
	phin1 += f;
	phin1.ExtrapGhost(0);
	F(phin1);
	phin1 += f;
	phin += phin1;
	phin *= 0.5;

	phin.ExtrapGhost(0);
}

void LevelSetSolver::TVDRK3()
{
	phin1 = phin;
	F(phin);
	phin1 += f;
	phin1.ExtrapGhost(0);
	F(phin1);
	phin1 += f;

	phin1 /= 3.0;
	phin1 += phin;
	phin1 *= 0.75;

	phin1.ExtrapGhost(0);
	F(phin1);
	phin1 += f;

	phin1 *= 2.0;
	phin += phin1;
	phin /= 3.0;

	phin.ExtrapGhost(0);
}

void LevelSetSolver::SimpleReinit()
{
	double ee = 9 * h_ * h_;
	double cfl = 0.5;
	double dt = cfl * h_, tmax = 0.8 * l_;

	// 这里用 f 替代 S
	for (int i = 0; i != n_; i++) {
		for (int j = 0; j != n_; j++) {
			f(i, j) = phin(i, j) / sqrt(phin(i, j) * phin(i, j) + ee);
		}
	}

	// pseudo time iteration
	double t = 0;
	double err = 1.0, errmin = 1e-4 * l_; // TODO: 收敛准则待考虑
	double a, b, c, d, ap, am, bp, bm, cp, cm, dp, dm, G;
	while (t<tmax && err>errmin)
	{
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {

				a = (phin(i, j) - phin(i - 1, j)) / h_;
				b = (phin(i + 1, j) - phin(i, j)) / h_;
				c = (phin(i, j) - phin(i, j - 1)) / h_;
				d = (phin(i, j + 1) - phin(i, j)) / h_;

				if (f(i, j) > 0) {
					ap = max(a, 0.0);
					bm = min(b, 0.0);
					cp = max(c, 0.0);
					dm = min(d, 0.0);

					G = sqrt(max(ap * ap, bm * bm) + max(cp * cp, dm * dm)) - 1.0;
				}
				else if (f(i, j) < 0)
				{
					am = min(a, 0.0);
					bp = max(b, 0.0);
					cm = min(c, 0.0);
					dp = max(d, 0.0);

					G = sqrt(max(am * am, bp * bp) + max(cm * cm, dp * dp)) - 1.0;
				}
				else {
					G = 0.0;
				}

				phin1(i, j) = phin(i, j) - dt * f(i, j) * G;
			}
		}

		// error calculation
		phin1 -= phin;
		err = 0.0;
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {
				err += phin1(i, j) * phin1(i, j);
			}
		}
		err = sqrt(err);
		// TEST
		//cout << err << endl;

		// update
		phin += phin1;
		phin.ExtrapGhost(0);

		t += dt;
	}
	std::cout << "Inner iteration steps: " << t / dt << endl;
}

void LevelSetSolver::SubcellFix()
{
	phin.ExtrapGhost();

	double cfl = 0.5;
	double dt = cfl * h_, tmax = 0.8 * l_;

	f = 1.0;
	D = 0.0;
	double phi0P, phi0W, phi0E, phi0S, phi0N;
	for (int i = 0; i != n_; i++) {
		for (int j = 0; j != n_; j++) {
			// sgn function
			if (phin(i, j) < 0) {
				f(i, j) = -1.0;
			}

			phi0P = phin(i, j);
			phi0W = phin(i - 1, j);
			phi0E = phin(i + 1, j);
			phi0S = phin(i, j - 1);
			phi0N = phin(i, j + 1);

			// Sigma
			if (phi0P * phi0W < 0 || phi0P * phi0E < 0 || \
				phi0P * phi0S < 0 || phi0P * phi0N < 0) {
				Sigma(i, j) = 1.0;
			}
			else {
				Sigma(i, j) = -1.0;
			}

			if (Sigma(i, j) > 0) {
				D(i, j) = max(abs(phi0E - phi0P), abs(phi0P - phi0W));
				D(i, j) = max(D(i, j), abs(phi0N - phi0P));
				D(i, j) = max(D(i, j), abs(phi0P - phi0S));
				D(i, j) = max(D(i, j), \
					0.5 * sqrt((phi0E - phi0W) * (phi0E - phi0W) + (phi0N - phi0S) * (phi0N - phi0S)));
				D(i, j) = h_ * phin(i, j) / D(i, j);
			}
		}
	}

	// pseudo time iteration
	double t = 0;
	double err = 1.0, errmin = 1e-10 * l_; // TODO: 收敛准则待考虑
	double a, b, c, d, ap, am, bp, bm, cp, cm, dp, dm, G;
	while (t<tmax && err>errmin)
	{
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {

				if (Sigma(i, j) > 0) {
					phin1(i, j) = phin(i, j) - dt / h_ * (f(i, j) * abs(phin(i, j)) - D(i, j));
				}
				else {
					a = (phin(i, j) - phin(i - 1, j)) / h_;
					b = (phin(i + 1, j) - phin(i, j)) / h_;
					c = (phin(i, j) - phin(i, j - 1)) / h_;
					d = (phin(i, j + 1) - phin(i, j)) / h_;

					if (f(i, j) > 0) {
						ap = max(a, 0.0);
						bm = min(b, 0.0);
						cp = max(c, 0.0);
						dm = min(d, 0.0);

						G = sqrt(max(ap * ap, bm * bm) + max(cp * cp, dm * dm)) - 1.0;
					}
					else if (f(i, j) < 0)
					{
						am = min(a, 0.0);
						bp = max(b, 0.0);
						cm = min(c, 0.0);
						dp = max(d, 0.0);

						G = sqrt(max(am * am, bp * bp) + max(cm * cm, dp * dp)) - 1.0;
					}
					else {
						G = 0.0;
					}

					phin1(i, j) = phin(i, j) - dt * f(i, j) * G;
				}
			}
		}

		// error calculation
		phin1 -= phin;
		err = 0.0;
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {
				err += phin1(i, j) * phin1(i, j);
			}
		}
		err = sqrt(err);
		// TEST
		//cout << err << endl;

		// update
		phin += phin1;
		phin.ExtrapGhost(0);

		t += dt;
	}
	std::cout << "Inner iteration steps: " << t / dt << endl;
}

void LevelSetSolver::Reinit()
{
	switch (scheme_r_)
	{
	case 1:
		SimpleReinit();
		break;
	case 2:
		SubcellFix();
		break;
	default: // no reinit
		break;
	}
}

LevelSetSolver::LevelSetSolver(int n, double dt, double tmax)
{
	n_ = n;
	h_ = 1.0 / double(n_);

	dt_ = dt;
	tmax_ = tmax;
	step_ = 0;
	stepmax_ = ceil(tmax / dt);

	phin.Init(n, n, 2, 0.0);
	phin1.Init(phin);
	f.Init(phin);
	u.Init(n, n, 0.0);
}

LevelSetSolver::LevelSetSolver(int n, double dt, double tmax, int scheme_s, int scheme_t, int scheme_r, int step_r)
{
	n_ = n;
	h_ = 1.0 / double(n_);

	dt_ = dt;
	tmax_ = tmax;
	step_ = 0;
	stepmax_ = ceil(tmax / dt);

	phin.Init(n, n, 2, 0.0);
	phin1.Init(phin);
	f.Init(phin);
	u.Init(n, n, 0.0);

	scheme_s_ = scheme_s;
	scheme_t_ = scheme_t;
	scheme_r_ = scheme_r;
	step_r_ = step_r;
}

LevelSetSolver::~LevelSetSolver()
{
}

void LevelSetSolver::Init(int n, double dt, double tmax, int scheme_s, int scheme_t, int scheme_r, int step_r)
{
	n_ = n;
	h_ = 1.0 / double(n_);

	dt_ = dt;
	tmax_ = tmax;
	step_ = 0;
	stepmax_ = ceil(tmax / dt);

	phin.Init(n, n, 2, 0.0);
	phin1.Init(phin);
	f.Init(phin);
	u.Init(n, n, 0.0);

	scheme_s_ = scheme_s;
	scheme_t_ = scheme_t;
	scheme_r_ = scheme_r;
	step_r_ = step_r;
}

void LevelSetSolver::InitVelocity(int type_init)
{
	double x, y;
	switch (type_init)
	{
	case 1: // streching
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {
				x = (i)*h_;
				y = (j + 0.5) * h_;
				u.w(i, j) = -2 * pi * cos(pi * (x - 0.5)) * sin(pi * (y - 0.5));

				x = (i + 1.0) * h_;
				y = (j + 0.5) * h_;
				u.e(i, j) = -2 * pi * cos(pi * (x - 0.5)) * sin(pi * (y - 0.5));

				x = (i + 0.5) * h_;
				y = (j)*h_;
				u.s(i, j) = -2 * pi * sin(pi * (x - 0.5)) * cos(pi * (y - 0.5));

				x = (i + 0.5) * h_;
				y = (j + 1.0) * h_;
				u.n(i, j) = -2 * pi * sin(pi * (x - 0.5)) * cos(pi * (y - 0.5));
			}
		}
		break;
	case 2: // tearing
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {
				x = (i)*h_;
				y = (j + 0.5) * h_;
				u.w(i, j) = sin(4.0 * pi * (x + 0.5)) * sin(4.0 * pi * (y + 0.5));

				x = (i + 1.0) * h_;
				y = (j + 0.5) * h_;
				u.e(i, j) = sin(4.0 * pi * (x + 0.5)) * sin(4.0 * pi * (y + 0.5));

				x = (i + 0.5) * h_;
				y = (j)*h_;
				u.s(i, j) = cos(4.0 * pi * (x + 0.5)) * cos(4.0 * pi * (y + 0.5));

				x = (i + 0.5) * h_;
				y = (j + 1.0) * h_;
				u.n(i, j) = cos(4.0 * pi * (x + 0.5)) * cos(4.0 * pi * (y + 0.5));
			}
		}
		break;
	case 3: // rotating
		for (int i = 0; i != n_; i++) {
			for (int j = 0; j != n_; j++) {
				y = (j + 0.5) * h_;
				u.w(i, j) = -2 * pi * (y - 0.5);
				u.e(i, j) = -2 * pi * (y - 0.5);

				x = (i + 0.5) * h_;
				u.s(i, j) = 2 * pi * (x - 0.5);
				u.n(i, j) = 2 * pi * (x - 0.5);
			}
		}
		break;
	default:
		break;
	}
}

void LevelSetSolver::InitPhi(int type_init)
{
	double x, y;
	double xc, yc, rc;

	Array in;

	Array Arc;
	Array LineSeg1;
	Array LineSeg2;
	Array Point1;
	Array Point2;
	Array Point3;
	double xp1;
	double yp1;
	double xp2;
	double yp2;
	double xp3;
	double yp3;

	double term = 0.0;

	switch (type_init)
	{
	case 1: // assignment 1 (Rider & Kothe variant problem)
		xc = 0.5;
		yc = 0.3;
		rc = 0.2;

		for (int i = -2; i != n_ + 2; i++) {
			for (int j = -2; j != n_ + 2; j++) {
				x = (i + 0.5) * h_;
				y = (j + 0.5) * h_;
				phin(i, j) = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)) - rc;
			}
		}
		break;


	case 2: // assignment 2 (Zalesak variant problem)
		xc = 0.5;
		yc = 0.5;
		rc = 0.4;

		in.Init(n_, n_, 2, 1.0);

		Arc.Init(n_, n_, 2, 10.0);
		LineSeg1.Init(n_, n_, 2, 10.0);
		LineSeg2.Init(n_, n_, 2, 10.0);
		Point1.Init(n_, n_, 2, 10.0);
		Point2.Init(n_, n_, 2, 10.0);
		Point3.Init(n_, n_, 2, 10.0);
		xp1 = xc - rc * cos(pi / 6.0);
		yp1 = yc + rc * sin(pi / 6.0);
		xp2 = xc;
		yp2 = yc;
		xp3 = xc - rc * cos(pi / 6.0);
		yp3 = yc - rc * sin(pi / 6.0);

		term = 0.0;
		for (int i = -2; i != n_ + 2; i++) {
			for (int j = -2; j != n_ + 2; j++) {
				x = (i + 0.5) * h_;
				y = (j + 0.5) * h_;

				Point1(i, j) = sqrt((x - xp1) * (x - xp1) + (y - yp1) * (y - yp1));
				Point2(i, j) = sqrt((x - xp2) * (x - xp2) + (y - yp2) * (y - yp2));
				Point3(i, j) = sqrt((x - xp3) * (x - xp3) + (y - yp3) * (y - yp3));

				term = y - x * tan(pi / 3.0);
				if (term >= yp2 - xp2 * tan(pi / 3.0) && term <= yp1 - xp1 * tan(pi / 3.0)) {
					term = y + x * tan(pi / 6.0);
					term -= yp2 + xp2 * tan(pi / 6.0);
					LineSeg1(i, j) = abs(term * cos(pi / 6.0));
				}

				term = y + x * tan(pi / 3.0);
				if (term <= yp2 + xp2 * tan(pi / 3.0) && term >= yp3 + xp3 * tan(pi / 3.0)) {
					term = y - x * tan(pi / 6.0);
					term -= yp2 - xp2 * tan(pi / 6.0);
					LineSeg2(i, j) = abs(term * cos(pi / 6.0));
				}

				term = (x - xc) / sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
				if (term > cos(pi * 5.0 / 6.0)) {
					term = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
					Arc(i, j) = abs(term - rc);
				}
			}
		}

		for (int i = -2; i != n_ + 2; i++) {
			for (int j = -2; j != n_ + 2; j++) {
				// distance function
				phin(i, j) = min(LineSeg1(i, j), LineSeg2(i, j));
				phin(i, j) = min(phin(i, j), Point1(i, j));
				phin(i, j) = min(phin(i, j), Point2(i, j));
				phin(i, j) = min(phin(i, j), Point3(i, j));
				phin(i, j) = min(phin(i, j), Arc(i, j));

				// sign distance function
				x = (i + 0.5) * h_;
				y = (j + 0.5) * h_;
				term = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
				if (term < rc) {

					if (y + x * tan(pi / 6.0) > yc + xc * tan(pi / 6.0) || \
						y - x * tan(pi / 6.0) < yc - xc * tan(pi / 6.0)) {
						phin(i, j) *= -1.0;
					}

				}
			}
		}

		break;


	case 3: // (Rider & Kothe 1995)
		xc = 0.5;
		yc = 0.75;
		rc = 0.15;

		for (int i = -2; i != n_ + 2; i++) {
			for (int j = -2; j != n_ + 2; j++) {
				x = (i + 0.5) * h_;
				y = (j + 0.5) * h_;
				phin(i, j) = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)) - rc;
			}
		}
		break;


	default: // assignment 1
		xc = 0.5;
		yc = 0.3;
		rc = 0.2;

		for (int i = -2; i != n_ + 2; i++) {
			for (int j = -2; j != n_ + 2; j++) {
				x = (i + 0.5) * h_;
				y = (j + 0.5) * h_;
				phin(i, j) = sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc)) - rc;
			}
		}
		break;
	}
}

void LevelSetSolver::Calculation()
{
	if (scheme_r_ == 2) { // for subcell fix
		Sigma.Init(n_, n_, 2, -1.0);
		D.Init(n_, n_, 2, 1e-10);
	}


	WriteAll();

	switch (scheme_t_)
	{
	case 1:
		while (step_ < stepmax_)
		{
			EULER();
			step_++;

			if (step_ % step_r_ == 0) {
				Reinit();
			}
		}
		break;
	case 2:
		while (step_ < stepmax_)
		{
			TVDRK2();
			step_++;

			if (step_ % step_r_ == 0) {
				Reinit();
			}
		}
		break;
	case 3:
		while (step_ < stepmax_)
		{
			TVDRK3();
			step_++;

			if (step_ % step_r_ == 0) {
				cout << "Step " << step_ << "\t t = " << step_ * dt_ << endl;
				Reinit();
				WriteAll();
			}
		}
		break;
	default:
		while (step_ < stepmax_)
		{
			TVDRK3();
			step_++;

			if (step_ % step_r_ == 0) {
				Reinit();
			}
		}
		break;
	}
	Reinit();
}

void LevelSetSolver::WriteVelocity()
{
	ofstream fout;
	string fname = "velocity_" + to_string(step_) + ".dat";
	fout.open(fname);

	// header
	fout << "VARIABLES = \"x\", \"y\", \"u\", \"v\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << n_ << endl;
	fout << "        J = " << n_ << endl;
	fout << "SOLUTIONTIME = " << step_ * dt_ << endl;

	// flow variables
	double x, y;

	for (int j = 0; j != n_; j++) {
		for (int i = 0; i != n_; i++) {

			x = (i + 0.5) * h_;
			y = (j + 0.5) * h_;

			fout << left << setw(8) << x << ' '							\
				<< left << setw(8) << y << ' '							\
				<< scientific << setprecision(12) << u.x(i, j) << ' '	\
				<< scientific << setprecision(12) << u.y(i, j) << endl;
		}
	}

	fout.close();
}

void LevelSetSolver::WriteVelocity(string fname)
{
	ofstream fout;
	fout.open(fname);

	// header
	fout << "VARIABLES = \"x\", \"y\", \"u\", \"v\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << n_ << endl;
	fout << "        J = " << n_ << endl;
	fout << "SOLUTIONTIME = " << step_ * dt_ << endl;

	// flow variables
	double x, y;

	for (int j = 0; j != n_; j++) {
		for (int i = 0; i != n_; i++) {

			x = (i + 0.5) * h_;
			y = (j + 0.5) * h_;

			fout << left << setw(8) << x << ' '							\
				<< left << setw(8) << y << ' '							\
				<< scientific << setprecision(12) << u.x(i, j) << ' '	\
				<< scientific << setprecision(12) << u.y(i, j) << endl;
		}
	}

	fout.close();
}

void LevelSetSolver::WritePhi()
{
	ofstream fout;
	string fname = "phi_" + to_string(step_) + ".dat";
	fout.open(fname);

	// header
	fout << "VARIABLES = \"x\", \"y\", \"phi\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << n_ << endl;
	fout << "        J = " << n_ << endl;
	fout << "SOLUTIONTIME = " << step_ * dt_ << endl;

	// flow variables
	double x, y;

	for (int j = 0; j != n_; j++) {
		for (int i = 0; i != n_; i++) {

			x = (i + 0.5) * h_;
			y = (j + 0.5) * h_;

			fout << left << setw(8) << x << ' '							\
				<< left << setw(8) << y << ' '							\
				<< scientific << setprecision(12) << phin(i, j) << endl;
		}
	}

	fout.close();
}

void LevelSetSolver::WritePhi(string fname)
{
	ofstream fout;
	fout.open(fname);

	// header
	fout << "VARIABLES = \"x\", \"y\", \"phi\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << n_ << endl;
	fout << "        J = " << n_ << endl;
	fout << "SOLUTIONTIME = " << step_ * dt_ << endl;

	// flow variables
	double x, y;

	for (int j = 0; j != n_; j++) {
		for (int i = 0; i != n_; i++) {

			x = (i + 0.5) * h_;
			y = (j + 0.5) * h_;

			fout << left << setw(8) << x << ' '							\
				<< left << setw(8) << y << ' '							\
				<< scientific << setprecision(12) << phin(i, j) << endl;
		}
	}

	fout.close();
}

void LevelSetSolver::WriteAll()
{
	ofstream fout;
	string fname = "output_" + to_string(step_) + ".dat";
	fout.open(fname);

	// header
	fout << "VARIABLES = \"x\", \"y\", \"u\", \"v\", \"phi\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << n_ << endl;
	fout << "        J = " << n_ << endl;
	fout << "SOLUTIONTIME = " << step_ * dt_ << endl;

	// flow variables
	double x, y;

	for (int j = 0; j != n_; j++) {
		for (int i = 0; i != n_; i++) {

			x = (i + 0.5) * h_;
			y = (j + 0.5) * h_;

			fout << left << setw(8) << x << ' '							\
				<< left << setw(8) << y << ' '							\
				<< scientific << setprecision(12) << u.x(i, j) << ' '	\
				<< scientific << setprecision(12) << u.y(i, j) << ' '	\
				<< scientific << setprecision(12) << phin(i, j) << endl;
		}
	}

	fout.close();
}

void LevelSetSolver::WriteAll(string fname)
{
	ofstream fout;
	fout.open(fname);

	// header
	fout << "VARIABLES = \"x\", \"y\", \"u\", \"v\", \"phi\"" << endl;
	fout << "ZONE    F = point" << endl;
	fout << "        I = " << n_ << endl;
	fout << "        J = " << n_ << endl;
	fout << "SOLUTIONTIME = " << step_ * dt_ << endl;

	// flow variables
	double x, y;

	for (int j = 0; j != n_; j++) {
		for (int i = 0; i != n_; i++) {

			x = (i + 0.5) * h_;
			y = (j + 0.5) * h_;

			fout << left << setw(8) << x << ' '							\
				<< left << setw(8) << y << ' '							\
				<< scientific << setprecision(12) << u.x(i, j) << ' '	\
				<< scientific << setprecision(12) << u.y(i, j) << ' '	\
				<< scientific << setprecision(12) << phin(i, j) << endl;
		}
	}

	fout.close();
}
