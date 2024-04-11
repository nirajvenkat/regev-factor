/// This Q# program sets a qubit register in a Gaussian superposition
/// using the arbitrary state preperation algorithm
/// from [Möttönen et al.](arXiv:quant-ph/0407010)
namespace Gaussian {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;

    @EntryPoint()
    operation GaussianStatePrepM() : Int {
        // Set number of qubits
        let n = 4;
        use reg = Qubit[n];
        let N = 2 ^ Length(reg);
        // Mean of the Gaussian
        let mu = N / 2;

        mutable amplitudes = [0.0, size = N];
        for i in IndexRange(amplitudes) {
            // 3σ corresponds to 99.7% of our kernel
            let sigma = (N - 1) / 3;
            let sigmaSq = IntAsDouble(sigma ^ 2);
            let prefactor = 1.0 / Sqrt(2.0 * PI() * sigmaSq);
            let weight = prefactor * E() ^ (-0.5 * IntAsDouble((i - mu) ^ 2) / sigmaSq);
            set amplitudes w/= i <- weight;
        }
        
        // Normalize PDF
        let ampNormed = PNormalized(2.0, amplitudes);
        MottonenStatePrep(ampNormed, reg);

        DumpMachine();
        
        // Measure the qubits.
        let m = MeasureInteger(reg);
        BoolArrayAsInt(Reversed(IntAsBoolArray(m, n)))
        // MeasureInteger(reg)
    }

    internal function GenAngles(x : Double[]) : Double[] {
        let N = Length(x);
        mutable angles = [0.0, size = N / 2];

        if N > 1 {
            mutable innerX = [0.0, size = N / 2];
            for k in IndexRange(innerX) {
                let s1 = AbsD(x[2 * k]) ^ 2.0;
                let s2 = AbsD(x[2 * k + 1]) ^ 2.0;
                set innerX w/= k <- Sqrt(s1 + s2);
            }

            let innerAngles = GenAngles(innerX);
            for k in IndexRange(innerX) {
                set angles w/= k <- 0.0;
                if innerX[k] != 0.0 {
                    let numer = x[2 * k + 1];
                    let denom = innerX[k];
                    let angleCandidate = 2.0 * ArcSin(numer / denom);
                    set angles w/= k <- (x[2 * k] > 0.0) ?
                                         angleCandidate |
                                         2.0 * PI() - angleCandidate;
                }
            }
            set angles = innerAngles + angles;
        }
        
        angles
    }

    operation MottonenStatePrep(amplitudes : Double[], reg : Qubit[]) : Unit is Adj + Ctl {
        let angles = GenAngles(amplitudes);
        for level in IndexRange(reg) {
            if level == 0 {
                Ry(angles[0], reg[level]);
            }
            else {
                let bitLen = (1 <<< level) - 1;
                for i in 0..bitLen {
                    let angleIdx = bitLen + i;
                    let bitArray = Reversed(IntAsBoolArray(i, level));
                    ApplyControlledOnBitString(bitArray, Ry, reg[0..level - 1], (angles[angleIdx], reg[level]));
                }
            }
        }
    }
}