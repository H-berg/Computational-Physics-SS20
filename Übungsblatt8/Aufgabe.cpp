#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>
#include "Eigen/Dense"
#define _USE_MATH_DEFINES

//using namespace std;
using namespace Eigen;
using namespace std;

//Initialisierung der Startparameter
//Geschwindigkeiten und Orte der Teilchen werden je in einer Matrix gespeichert, in welcher jede Spalte für ein Teilchen steht
void Initialisierung(MatrixXd &r, MatrixXd &v, const int N, const double L, double T)
{
    int count = 0;

    //Ortsvektor r
    for (int n = 0; n < sqrt(N); n++)
    {
        for (int m = 0; m < sqrt(N); m++)
        {

            r(0, count) = (1 + 2 * n) * L / 8;
            r(1, count) = (1 + 2 * m) * L / 8;
            count++;
        }
    }

    //Es werden zufällige Geschwindigkeiten gewählt
    //und der Mittelwerte dieser abgezogen um die Schwerpunktsgeschwindigkeit auf 0 zu setzen
    v = MatrixXd::Random(2, N);
    Vector2d mu = v.rowwise().sum() / N;
    for (int i = 0; i < N; i++)
    {
        v.col(i) = (v.col(i) - mu);
    }


    //Umskalierung der Geschwindigkeiten
    double T2 = v.colwise().squaredNorm().sum() / (2 * N - 2); // m = 1 k_B = 1
    double n = sqrt(T / T2);
    v = v * n;

}

// gibt den Vektor zum nächstgelegenden Teilchen wieder (Bild- oder Realteilchen)
// Ist weder Bildteilchen oder Realteilchen innerhalb des Radius L/2 (eines Teilchens) wird 0 als Abstand zurückgegeben
// Die Kraft und Potentialfunktion sind wiederum so programmiert, dass diese für r=0 ebenfalls 0 zurückgeben, da es dann keine WW zwischen den Teilchen gibt
VectorXd naechstesTeilchen(VectorXd r1, VectorXd r2, const double L)
{
    VectorXd Ln(2);
    //Es wird Realteilchen und jedes Bildteilchen überprüft
    //höchstens ein Teilchen kann sich im Radius L/2 befinden
    //wird dieses gefunden wird die Funktion durch return abgebrochen
    for (int i = -1; i <= 1; i++)
    {
        for (int j = -1; j <= 1; j++)
        {
            Ln(0) = i * L;
            Ln(1) = j * L;
            VectorXd abstand = r1 - (r2 + Ln);
            if (abstand.norm() < L / 2)
            {
                return abstand;
            }
        }
    }
    return VectorXd::Zero(2);
}

//Kraft des LJ-Potentials
VectorXd LJ_Kraft(VectorXd r, const double L)
{
    if (r.norm()< 1E-3)
    {
        return VectorXd::Zero(2);
    }
    VectorXd ret = r * 24 * (2 * pow(r.norm(), -14) - pow(r.norm(), -8));
    return ret;
}

//LJ-Potential
double LJ_Pot(VectorXd r, const double L)
{
    double rn = r.norm();
    if(rn < 1E-3)
    {
        return 0;
    }
    rn = pow(rn, -6); //So reduziert sich die Berechnung von r^12 auf eine quadrierung
    double ret = 4 * (pow(rn, 2) - rn) - 4 * (pow(L / 2, -12) - pow(L / 2, -6)); // Potential soll stetig sein somit wird es um r_c=L/2 verschoben
    return ret;
}

//Berechnet die Kraftmatrix mit der Kraft auf das entsprechende Teilchen in der entsprechenden Spalte
MatrixXd Kraft(MatrixXd &r, const double L, const int N)
{
    MatrixXd Kraft = MatrixXd::Zero(2, N);
    VectorXd r_next(2), F_temp(2);

    for (int i = 0; i < (N - 1); i++) // es wird die Kraft auf Teilchen i berechnet, die Kraft des letzten Teilchen wird bereits mit den vorherigen Rechnungen komplett berechnet
    {

        for (int n = (i + 1); n < N; n++) //indem die Wechselwirkung mit den anderen Teilchen addiert wird, die Kraft für die Teilchen j < i wurde bereits addiert;
        {
            r_next = naechstesTeilchen(r.col(i), r.col(n), L); // Abstand zum nächsten Teilchen j
            F_temp = LJ_Kraft(r_next, L);                      //Kraft zwischen den beiden Teilchen
            Kraft.col(i) += F_temp;
            Kraft.col(n) -= F_temp; // die Kraft von i auf j ist die negative Kraft von j auf i, so spart man die Hälfte der Rechnungen
        }
    }

    return Kraft;
}

//Die Funktion paast den Ortsvektor entsprechend den Randbedingungen an
void Randbedingungen(MatrixXd &r, const double L, const int N)
{
    VectorXd d = VectorXd::Zero(2);

//Vermutlich kommt es hier bei T=100 zu Problemen weil die Geschwindigkeit gegen Unendlich geht u.Ä.
    for (int i = 0; i < N; i++)
    {
        d(0) = floor(r(0, i) / L);
        d(1) = floor(r(1, i) / L); //Verwendung der Floorfunktion wie in der Vorlesung

        r.col(i) -= (L * d);
    }
}

//Implementierung des Verlet Algorithmus
void Verlet_Iteration(MatrixXd &r, MatrixXd &v, MatrixXd &a, const double h, const double L, const int N)
{
    MatrixXd a_temp = a; // m = 1

    r = r + v * h + a * h * h / 2;

    //Berücksichtigung der periodischen Randbedingungen
    Randbedingungen(r, L, N);

    a = Kraft(r, L, N);
    v = v + (a + a_temp) * h / 2;
}

//Schwerpunktsgeschwindigkeit
VectorXd SP_Geschw(MatrixXd v, const int N)
{
    VectorXd ret;
    ret = v.colwise().sum() / N;
    return ret;
}

//Temperatur
double Temperatur(const int N, double E_kin)
{
    double ret;
    ret = 2 * E_kin / (2 * N - 2);
    return ret;
}

//potentielle Energie
double pot_Energie(MatrixXd r, const double N, const double L)
{
    VectorXd r0;
    double pot = 0;

    for (int i = 0; i < N - 1; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            r0 = naechstesTeilchen(r.col(i), r.col(j), L);
            pot += LJ_Pot(r0, L);
        }
    }

    return pot;
}

//Gibt die Anzahl der Teilchen pro Bin zurück
VectorXi TeilchenInBin(MatrixXd r, const int N_h, const double L, const int N)
{
    double Delta_r = L / (2 * N_h); //Dicke des Kreisrings

    VectorXd Ln(2), r1(2), r2(2);
    VectorXi bins = VectorXi::Zero(N_h);

    for (int l = 0; l < N_h; l++) //es wird Bin l betrachtet
    {
        double r_min = Delta_r * l;
        double r_max = Delta_r * (l + 1);
        for (int m = 0; m < (N-1); m++)
        {
            r1 = r.col(m);
            for (int n = (m+1); n < N; n++)//ähnlich wie bei der Kraft wird versucht, die gleichen Betrachtungen nicht doppelt durchzuführen
            {
                r2 = r.col(n);
                //Es werden Teilchen und Bildteilchen überprüft
                for (int i = -1; i <= 1; i++)
                {
                    for (int j = -1; j <= 1; j++) 
                    {
                        Ln(0) = i * L;
                        Ln(1) = j * L;
                        VectorXd abstand = r1 - (r2 + Ln);
                        if (abstand.norm() < r_max && abstand.norm() >= r_min)
                        {
                            bins(l)+=2; // + zwei da wenn n im l-ten Bin von m dann ist auch m in l-ten Bin von n
                            i=1;
                            break;
                            //Sobald von einem Teilchen, Real/ oder Bildteilchen im Bin ist kann die Iteration abgebrochen werden
                            //da höchstens für ein Teilchen r<L/2 gilt. Somit spart man unnötige Iterationsschritte
                        }
                    }
                }
            }
        }
    }
    return bins;
}

//kinetische Energie
double kin_Energie(MatrixXd v, const int N)
{
    double r;
    r = v.colwise().squaredNorm().sum() / 2;
    return r;
}

//Simuliert das Verhalten bis der Gleichgewichtszustand eintritt
void Aequilibrierung(MatrixXd &r, MatrixXd &v, const int N, const double L, double tmax, ofstream &f, string Name)
{
    double t = 0, h = 0.01;
    MatrixXd a = Kraft(r, L, N);
    VectorXd SP_G;
    double T, E_pot, E_kin;

    f.open(Name + ".txt");
    f << "# t \t SP-Geschwindigkeit x, y \t E_kin \t E_pot \t Temperatur" << endl;

    //Verlet Algorithmus, gleichzeitig werden auch die entsprechenden Größen (Aufgabeteil b) berechnet und gespeichert
    while (t < tmax)
    {
        t += h;
        Verlet_Iteration(r, v, a, h, L, N); // update von a, r und v
        SP_G = SP_Geschw(v, N);
        E_pot = pot_Energie(r, N, L);
        E_kin = kin_Energie(v, N);
        T = Temperatur(N, E_kin);

        f << t << "\t" << SP_G(0) << "\t" << SP_G(1) << "\t" << E_kin << "\t" << E_pot << "\t" << T << "\n";
    }
    f.flush();
    f.close();
}

//Berechnung der Paarkorrelationsfunktion
VectorXd Paarkorrelationsfunktion(VectorXi P, const int N_h, double L , const int N, int n)
{
    VectorXd g = VectorXd::Zero(N_h);
    P /= (n-2); //Zeitmittellung durch Teilen durch die Anzahl an Zeitschritten
    double rho = N/(L*L), Delta_V;
    double Delta_r = L / (2 * N_h);

    for(int l = 0; l < N_h; l++)
    {
        //Für jeden Bin muss V berechnet werden
        Delta_V = M_PI*((Delta_r*(l+1))*(Delta_r*(l+1))- (Delta_r*l)*(Delta_r*l)); // Fläche des Kreisrings
        g(l)= P(l)/(N*rho*Delta_V);
    }
    return g;
}


void Messung(MatrixXd &r, MatrixXd &v, const int N, const int N_h, const double L, int n, ofstream &f, string Name) // n ist die # an Zeitschritten
{
    double t = 0, h = 0.01, E_kin;
    double tmax = n*h;
    MatrixXd a = Kraft(r, L, N);
    VectorXd T = VectorXd::Zero(n);
    unsigned int count = 0;
    VectorXi P = VectorXi::Zero(N_h);

    //Verlet Algorithmus
    while (t < tmax)
    {
        t += h;
        Verlet_Iteration(r, v, a, h, L, N); // update von a, r und v
        E_kin = kin_Energie(v, N);
        T[count] = Temperatur(N, E_kin);
        count ++;

        // Für die Paarkorrelationsfunktion wird die Anzahl der Teilchen pro Bin benötigt
        P += TeilchenInBin(r, N_h, L, N);

        //f << t << "\t" << SP_G(0) << "\t" << SP_G(1) << "\t" << E_kin << "\t" << E_pot << "\t" << T << "\n";
    }

    VectorXd g = Paarkorrelationsfunktion(P, N_h, L, N, n);
    double Delta_r = L / (2 * N_h);

    f.open(Name+".txt");
    double T_mean = T.mean();
    f << "# Starttemperatur " + Name << endl;
    f << "#Endtemperatur gemittelt über " << n << " Zeitschritte: " << T_mean << "\n\n";
    f << "#Radius \t Paarkorrelation \n";
    for(int i = 0; i < N_h; i++)
    {
        f << Delta_r*i << "\t" << g(i) <<"\n";
    }
    f.flush();
    f.close();
}

int main()
{
    //Initialisierung
    const int N = 16, n = pow(10, 4), N_h = 50; 
    const double L = 8, t_eq = 20; // Kastenmaß, Zeitschritte bis zur Äquilibrierung, Anzahl der Zeitschritte bei Messung
    double T = 1;
    MatrixXd r = MatrixXd::Zero(2, N);
    MatrixXd v = MatrixXd::Zero(2, N);

    //Aufgabenteil b)
    Initialisierung(r, v, N, L, T);
    ofstream f;
    Aequilibrierung(r, v, N, L, t_eq, f, "Aequilibrierung");

    //Aufgabenteil c)


    string TstringVec[3] = {"0.01", "1", "100"}; // 100 wird erstmal ausgelassen, da es hier zu numerischen Problemen kommt
    for (auto &Tstring : TstringVec)
    {
        r = MatrixXd::Zero(2, N);
        v = MatrixXd::Zero(2, N);
        const double T_start = stod(Tstring);
        Initialisierung(r, v, N, L, T_start);
        cout << T_start << endl;
        ofstream g;
        Aequilibrierung(r, v, N, L, t_eq, g, "Aequilibrierung"+Tstring);
        Messung(r, v, N, N_h, L, n, g, "Messung"+Tstring);
    }

    return 0;
}