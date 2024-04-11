/// This is a resource estimation Q# program for the quantum part of Regev's algorithm
/// that uses Fibonacci representation in Zeckendorf form for exponentiation.
namespace Regev {
    open Microsoft.Quantum.ResourceEstimation;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Random;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Unstable.Arithmetic;
    open Microsoft.Quantum.Unstable.StatePreparation;

    /// Choice of superposition
    /// 1: Gaussian
    /// 2: Uniform
    internal function SuperpositionState_() : Int { 1 }

    @EntryPoint()
    operation EstimateFibonacciExponentiation() : Int[] {
        // let N = 143; // 11*13
        let N = 8051; // 83*97
        
        // Constants in Regev paper
        let n = BitSizeI(N);
        let d = Floor(Sqrt(IntAsDouble(n)));
        let C = 2;
        let R = 2 ^ (C * d);
        // let D = 2 ^ Ceiling(Lg(2.0 * Sqrt(IntAsDouble(d)) * IntAsDouble(R))); // Ekera Dlog paper
        // let D = 2 ^ ((C + 2) * d); // RV optimized paper
        let D = 20; // Artificial clamped value

        // Choose superposition
        use ctl = Qubit[D];
        if SuperpositionState_() == 1 {
            GaussianStatePrep(ctl);
        } else {
            ApplyToEach(H, ctl);
        }
        
        mutable x = d ^ 2;
        while ApproxPrimeCounting(x) < d {
            set x = x <<< 1;
        }
        let a = Mapped(i -> i ^ 2, SieveEratosthenes(x));
        Fact(Length(a) >= d, "Error: need atleast d prime squares to continue");

        let zRange = 2 ^ ((n / d) + d);
        use tgt = Qubit[n];
        mutable res : Int[] = [];
        
        within {
            // In Regev's algorithm we sample the dual lattice vector d + 4 times.
            // Instruct resource estimator to simulate this repeated sampling.
            RepeatEstimates(d + 4);
        }
        apply {
            // Initialize product register to 1
            X(tgt[0]);

            // Exponentiate the product register
            for i in 0..d - 1 {
                let z_i = DrawRandomInt(0, zRange);
                FibonacciExponentiation(ctl, a[i], z_i, N, tgt);
            }
            
            Adjoint ApplyQFT(ctl);

            set res += [MeasureInteger(ctl)];
            ResetAll(tgt);
        }

        // Returns noisy samples from the dual lattice L*
        res
    }

    /// # Summary
    /// Prepare a discrete Gaussian superposition.
    operation GaussianStatePrep(reg : Qubit[]) : Unit {
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
    }

    /// # Summary
    /// Approximation of prime counting function 𝝅(x).
    ///
    /// # Description
    /// Lower bound on the number of primes less than or equal to `x`.
    function ApproxPrimeCounting(x: Int) : Int {
        let xd = IntAsDouble(x);
        let lxd = Log(xd);

        Round(xd / lxd)
    }

    /// # Summary
    /// Simple sieve to extract first `N` primes in time O(N log₂log₂N).
    function SieveEratosthenes(N : Int) : Int[] {
        mutable sieve = [true, size = N + 1];
        set sieve w/= 0 <- false;
        set sieve w/= 1 <- false;
        mutable i = 2;
        while i * i <= N {
            if sieve[i] {
                mutable k = i ^ 2;
                while k <= N {
                    set sieve w/= k <- false;
                    set k += i;
                }
            }
            set i += 1;
        }
        let sieveMapped = MappedByIndex((index, element) -> element ? index | -1, sieve);
        let sieveFiltered = Filtered(x -> x != -1, sieveMapped);

        sieveFiltered
    }
    
    /// # Summary
    /// Computes an exponential by scanning a Fibonacci representation in “low-to-high” order.
    /// This is the algorithm FibExp from [Kaliski](arXiv:1711.02491).
    /// See also Knuth: Appl. Math. Lett. Vol. I, No. 2, pp. III-VI, 1988
    ///
    /// # Note
    /// This is a classically assisted technique, not a full quantum-quantum
    /// multiplication which may be required in practice.
    operation FibonacciExponentiation(
        ctl : Qubit[],
        a : Int,
        k : Int,
        modulus : Int,
        out : Qubit[]
    ) : Unit is Adj + Ctl {
        let fRep = FibonacciRepresentation(k);
        for i in IndexRange(fRep) {
            if fRep[i] {
                Controlled ModularMultiplyByConstant(
                    ctl,
                    (
                        modulus,
                        ExpModI(a, FibonacciNumber(i + 2), modulus),
                        out
                    )
                )
            }
        }
    }

    /// # Summary
    /// Returns Fibonacci representation of `N` in time O(log₂log₂N).
    function FibonacciRepresentation(N : Int) : Bool[] {
        Fact(N > 0, "Number must be positive.");
        if N == 1 {
            return [true];
        }
        
        // One-sided binary search
        mutable r = 1;
        while FibonacciNumber(r) < N {
            set r <<<= 1;
        }

        // Once we have the upper bound, run binary search that
        // terminates at the largest Fibonacci number <= Fib(N).
        mutable l = r >>> 1;
        mutable found = false;
        while l < r and not found {
            let m = l + ((r - l) >>> 1);
            let f1 = FibonacciNumber(m);
            let f2 = FibonacciNumber(m + 1);

            if (f1 < N and f2 > N) or f1 == N {
                set l = m;
                set found = true;
            }
            elif f1 < N {
                set l = m + 1;
            }
            else {
                set r = m - 1;
            }
        }

        mutable repr = [false, size = l + 1];
        mutable n = N;
        while n > 0 {
            set n -= FibonacciNumber(l);
            set repr w/= l <- true;
            while l >= 1 and FibonacciNumber(l) > n {
                set l -= 1;
            }
        }

        repr[2 ...]
    }

    /// # Summary
    /// Computes the Nth Fibonacci number in O(log₂N) time
    /// using the [matrix method](https://wp.me/pZjY-1U).
    ///
    /// # Note
    /// [[ 1 1 ] ^n    [[ a b ]     [[ Fib(n+1)  Fib(n)   ]
    ///  [ 1 0 ]]   =   [ b c ]] =   [ Fib(n)    Fib(n-1) ]]
    /// See Eq. 8 [Dasdan](arXiv:1803.07199) for the recursive formula
    /// which is rearranged for the iterative implementation.
    function FibonacciNumber(N: Int) : Int {
        Fact(N >= 0, "Number must be non-negative.");
        mutable (a, b, c) = (1, 0, 1);
        for bit in Reversed(IntAsBoolArray(N, BitSizeI(N))) {
            if bit {
                set (a, b) = ((a + c) * b, b * b + c * c);
            } else {
                set (a, b) = (a * a + b * b, (a + c) * b);
            }
            set c = a + b;
        }
        
        b
    }
    
    /// # Summary
    /// Performs modular in-place multiplication by a classical constant.
    ///
    /// # Description
    /// Given the classical constants `c` and `modulus`, and an input quantum
    /// register |𝑦⟩ in little-endian format, this operation computes
    /// `(c*x) % modulus` into |𝑦⟩.
    ///
    /// # Input
    /// ## modulus
    /// Modulus to use for modular multiplication
    /// ## c
    /// Constant by which to multiply |𝑦⟩
    /// ## y
    /// Quantum register of target
    internal operation ModularMultiplyByConstant(modulus : Int, c : Int, y : Qubit[])
    : Unit is Adj + Ctl {
        use qs = Qubit[Length(y)];
        for idx in IndexRange(y) {
            let shiftedC = (c <<< idx) % modulus;
            Controlled ModularAddConstant(
                [y[idx]],
                (modulus, shiftedC, qs));
        }
        for idx in IndexRange(y) {
            SWAP(y[idx], qs[idx]);
        }
        let invC = InverseModI(c, modulus);
        for idx in IndexRange(y) {
            let shiftedC = (invC <<< idx) % modulus;
            Controlled ModularAddConstant(
                [y[idx]],
                (modulus, modulus - shiftedC, qs));
        }
    }

    /// # Summary
    /// Performs modular in-place addition of a classical constant into a
    /// quantum register.
    ///
    /// Given the classical constants `c` and `modulus`, and an input quantum
    /// register |𝑦⟩ in little-endian format, this operation computes
    /// `(x+c) % modulus` into |𝑦⟩.
    ///
    /// # Input
    /// ## modulus
    /// Modulus to use for modular addition
    /// ## c
    /// Constant to add to |𝑦⟩
    /// ## y
    /// Quantum register of target
    internal operation ModularAddConstant(modulus : Int, c : Int, y : Qubit[])
    : Unit is Adj + Ctl {
        body (...) {
            Controlled ModularAddConstant([], (modulus, c, y));
        }
        controlled (ctrls, ...) {
            // We apply a custom strategy to control this operation instead of
            // letting the compiler create the controlled variant for us in
            // which the `Controlled` functor would be distributed over each
            // operation in the body.
            //
            // Here we can use some scratch memory to save ensure that at most
            // one control qubit is used for costly operations such as
            // `AddConstant` and `CompareGreaterThenOrEqualConstant`.
            if Length(ctrls) >= 2 {
                use control = Qubit();
                within {
                    Controlled X(ctrls, control);
                } apply {
                    Controlled ModularAddConstant([control], (modulus, c, y));
                }
            } else {
                use carry = Qubit();
                Controlled IncByI(ctrls, (c, y + [carry]));
                Controlled Adjoint IncByI(ctrls, (modulus, y + [carry]));
                Controlled IncByI([carry], (modulus, y));
                Controlled ApplyIfLessOrEqualL(ctrls, (X, IntAsBigInt(c), y, carry));
            }
        }
    }
}
