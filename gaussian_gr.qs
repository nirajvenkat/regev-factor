/// This Q# program sets a qubit register in a Gaussian superposition
/// using the [Grover-Rudolph algorithm](arXiv:quant-ph/0208112)
/// implemented by [Klco et. al.](arXiv:1904.10440)
namespace Gaussian {
    open Microsoft.Quantum.ResourceEstimation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;

    @EntryPoint()
    operation GaussianStatePrepGR() : Int {
        // Set number of qubits
        let n = 4;
        use reg = Qubit[n];
        
        GroverRudolphStatePrep(reg);

        DumpMachine();
        
        // Measure the qubits.
        let m = MeasureInteger(reg);
        BoolArrayAsInt(Reversed(IntAsBoolArray(m, n)))
    }

    operation GroverRudolphStatePrep(reg : Qubit[]) : Unit is Adj + Ctl {
        let angles = GenAngles(Length(reg));
        Ry(angles[0][0], reg[0]);
        for l in 1..Length(reg) - 1 {
            for i in 0..(1 <<< l) - 1 {
                let bitArray = Reversed(IntAsBoolArray(i, l));
                ApplyControlledOnBitString(bitArray, Ry, reg[0..l-1], (angles[l][i], reg[l]));
            }
        }
    }
    
    internal function GenAngles(n : Int) : Double[][] {
        mutable angles : Double[][] = [[], size = n];
        for l in 0..n - 1 {
            mutable anglesInner : Double[] = [];
            for k in 0..2 ^ l - 1 {
                set anglesInner += [Theta(n, l, k)];
            }
            set angles w/= l <- anglesInner;
        }

        angles
    }

    internal function Theta(n : Int, l : Int, k : Int) : Double {
        let xMin = (1 <<< (n - l - 1)) * (2 * k + 1);
        let xMax = (1 <<< (n - l)) * (k + 1) - 1;
        let yMin = k * (1 <<< (n - l));
        let yMax = xMin - 1;
        
        mutable ns : Double[] = [];
        for x in xMin..xMax {
            set ns += [Psi(n, x) ^ 2.0];
        }
        mutable ds : Double[] = [];
        for y in yMin..yMax {
            set ds += [Psi(n, y) ^ 2.0];
        }

        let sumN = Fold((x, y) -> x + y, 0.0, ns);
        let sumD = Fold((x, y) -> x + y, 0.0, ds);

        2.0 * ArcTan(Sqrt(sumN / sumD))
    }

    internal function Psi(n : Int, x : Int) : Double {
        let N = 2 ^ n;
        let muC = 1.0 - (1.0 / IntAsDouble(N));
        let mu = IntAsDouble(2 ^ (n - 1)) * muC;
        let sigmaC = 0.5;
        let sigma = IntAsDouble(2 ^ (n - 1)) * sigmaC;
        let sigmaSq = sigma ^ 2.0;

        let xTerm = E() ^ (-((IntAsDouble(x) - mu) ^ 2.0) / (2.0 * sigmaSq));
        mutable allTerms : Double[] = [];
        for i in 0..N - 1 {
            set allTerms += [E() ^ (-((IntAsDouble(i) - mu) ^ 2.0) / (2.0 * sigmaSq))];
        }
        let norm = Sqrt(Fold((x, y) -> x + y, 0.0, allTerms));

        xTerm / norm
    }
}