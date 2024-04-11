/// This Q# program sets a qubit register in a Gaussian superposition
namespace Gaussian {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Unstable.StatePreparation;

    @EntryPoint()
    operation GaussianStatePrep() : Int {
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
        PreparePureStateD(ampNormed, reg);

        DumpMachine();
        
        // Measure the qubits.
        let m = MeasureInteger(reg);
        BoolArrayAsInt(Reversed(IntAsBoolArray(m, n)))
    }
}