# rustplex

A fast persistent homology library written in Rust.

---

## What is Persistent Homology?

Persistent homology is a method from **Topological Data Analysis (TDA)** that measures the "shape" of data across multiple scales. The core idea: as you grow a parameter ε (think of it as a proximity radius), topological features — connected components, loops, voids — appear and disappear. Persistence tracks *when* each feature is born and when it dies.

---

## Mathematical Foundations

### 1. Simplicial Complexes

A **simplicial complex** K is a combinatorial space built from simplices:

- A **0-simplex** is a vertex: {v}
- A **1-simplex** is an edge: {v₀, v₁}
- A **2-simplex** is a triangle: {v₀, v₁, v₂}
- A **k-simplex** is a filled k-dimensional polytope: {v₀, …, vₖ}

The complex must be closed under faces: if σ ∈ K and τ ⊆ σ, then τ ∈ K.

**Example.** Given four points {a, b, c, d}, one valid simplicial complex is:

```
K = { {a}, {b}, {c}, {d},
      {a,b}, {b,c}, {a,c},
      {a,b,c} }
```

This represents a filled triangle with an isolated vertex d.

---

### 2. The Vietoris–Rips Complex

Given a point cloud X ⊂ ℝⁿ and a scale parameter ε > 0, the **Vietoris–Rips complex** Rips(X, ε) contains a simplex {x₀, …, xₖ} whenever every pair of points in that set is within distance ε of each other:

```
{x₀, …, xₖ} ∈ Rips(X, ε)  ⟺  d(xᵢ, xⱼ) ≤ ε  for all i, j
```

**Example.** Five points arranged roughly in a circle. At small ε, you have five isolated vertices. As ε grows, edges appear between nearby points. At some ε*, a loop forms (a 1-cycle). At larger ε, triangles fill in and the loop is killed. The loop's birth and death times record its persistence.

---

### 3. Chain Groups and Boundary Operators

For each dimension k, the **k-th chain group** Cₖ(K; 𝔽) is the vector space over a field 𝔽 (commonly 𝔽₂ = {0,1}) spanned by the k-simplices of K. Elements are formal sums of simplices (chains).

The **boundary operator** ∂ₖ : Cₖ → Cₖ₋₁ maps each simplex to its oriented boundary:

```
∂ₖ({v₀, v₁, …, vₖ}) = Σᵢ (-1)ⁱ {v₀, …, v̂ᵢ, …, vₖ}
```

where v̂ᵢ means the i-th vertex is omitted.

**Example.** The boundary of the triangle {a, b, c} is:

```
∂₂({a,b,c}) = {b,c} − {a,c} + {a,b}
```

Over 𝔽₂ (where −1 = 1):

```
∂₂({a,b,c}) = {b,c} + {a,c} + {a,b}
```

The boundary of this boundary is zero: ∂₁ ∘ ∂₂ = 0. This is always true: **∂² = 0**.

---

### 4. Homology Groups

The **k-th homology group** measures k-dimensional holes:

```
Hₖ(K; 𝔽) = ker(∂ₖ) / im(∂ₖ₊₁)
```

- **ker(∂ₖ)** = k-cycles: chains with no boundary (closed loops, closed surfaces, …)
- **im(∂ₖ₊₁)** = k-boundaries: chains that bound something one dimension higher

A cycle that is *not* a boundary represents a genuine hole.

**Example.** For a circle (triangulated as three edges forming a loop):

```
H₀ = 𝔽    (one connected component)
H₁ = 𝔽    (one loop)
H₂ = 0    (no enclosed volume)
```

For a filled triangle:

```
H₀ = 𝔽    (one connected component)
H₁ = 0    (the loop is filled — it bounds the triangle)
```

The **Betti numbers** βₖ = dim(Hₖ) count the number of independent k-dimensional holes.

---

### 5. Filtrations

A **filtration** of K is a nested sequence of subcomplexes indexed by a parameter:

```
∅ = K₀ ⊆ K₁ ⊆ K₂ ⊆ … ⊆ Kₙ = K
```

For the Vietoris–Rips complex, the filtration is induced by the pairwise distances: each simplex enters at the scale equal to its diameter (the maximum pairwise distance among its vertices).

**Example.** Four points with pairwise distances:

```
d(a,b) = 1,  d(b,c) = 1,  d(c,d) = 1,  d(a,d) = 1
d(a,c) = √2, d(b,d) = √2
```

Filtration:
- ε < 1: four isolated vertices
- ε = 1: edges {a,b}, {b,c}, {c,d}, {a,d} appear — a square loop is born
- ε = √2: diagonals {a,c}, {b,d} appear, then triangles — the loop is killed

The loop lives from ε = 1 to ε = √2. Its persistence is √2 − 1.

---

### 6. Persistent Homology and the Reduction Algorithm

Persistent homology tracks how homology classes are born and die as the filtration parameter grows. The output is a set of **persistence pairs** (b, d): a class born at scale b and killed at scale d.

The key computational result (Zomorodian & Carlsson, 2005) is that persistence pairs can be read off from the **reduced boundary matrix**.

Given the filtration order σ₁, σ₂, …, σₙ on simplices, form the boundary matrix M where column j encodes ∂(σⱼ) as a binary vector over the filtration-ordered simplices. Then apply column reduction:

```
For j from 1 to n:
    while low(j) = low(i) for some i < j:
        add column i to column j (over 𝔽₂)
```

where **low(j)** is the index of the lowest nonzero entry in column j. After reduction:

- A zero column j: σⱼ creates a new homology class (it is a cycle not killed yet)
- A nonzero column j with low(j) = i: the pair (i, j) records birth at σᵢ and death at σⱼ

**Example.** Triangle {a, b, c} with edges and vertices in filtration order:

```
σ₁={a}, σ₂={b}, σ₃={c}, σ₄={a,b}, σ₅={b,c}, σ₆={a,c}, σ₇={a,b,c}
```

Boundary matrix (columns = simplices, rows = lower-dim simplices, 𝔽₂):

```
        σ₁  σ₂  σ₃  σ₄  σ₅  σ₆  σ₇
σ₁  [   0   0   0   1   0   1   0  ]
σ₂  [   0   0   0   1   1   0   0  ]
σ₃  [   0   0   0   0   1   1   0  ]
σ₄  [   0   0   0   0   0   0   1  ]
σ₅  [   0   0   0   0   0   0   1  ]
σ₆  [   0   0   0   0   0   0   1  ]
```

After reduction, the algorithm pairs (σ₁,σ₄), (σ₂,σ₅), and (σ₃,σ₆) — the three edges kill the three extra connected components. Column σ₇ reduces to zero after killing the 1-cycle {a,b}+{b,c}+{a,c}, so the loop (σ₆) is paired with (σ₇). The result is exactly the expected homology: one persistent H₀ class and one H₁ class born and immediately killed.

---

### 7. Persistence Diagrams and Stability

The output of the algorithm is a **persistence diagram** Dgm(K): a multiset of points (b, d) in the extended half-plane b ≤ d. Each point represents a topological feature. Features far from the diagonal (large d − b) are "significant"; those near the diagonal are noise.

**Bottleneck distance** between two diagrams D₁ and D₂:

```
d_B(D₁, D₂) = inf_{η: D₁→D₂} sup_{x ∈ D₁} ‖x − η(x)‖_∞
```

where the infimum is over all matchings η between the diagrams (with diagonal points included to allow unmatched features).

**Stability Theorem** (Cohen-Steiner, Edelsbrunner & Harer, 2007): For tame functions f, g on a triangulable space:

```
d_B(Dgm(f), Dgm(g)) ≤ ‖f − g‖_∞
```

Small perturbations in input data produce small changes in persistence diagrams. This is what makes persistent homology a robust descriptor for data analysis.

---

## Key References

- Edelsbrunner & Harer, *Computational Topology* (2010)
- Zomorodian, *Topology for Computing* (2005)
- Zomorodian & Carlsson, "Computing Persistent Homology," *Discrete & Computational Geometry* (2005)
- Edelsbrunner, Letscher & Zomorodian, "Topological Persistence and Simplification," *FOCS* (2002)
- Cohen-Steiner, Edelsbrunner & Harer, "Stability of Persistence Diagrams," *DCG* (2007)
- Bauer, "Ripser: Efficient Computation of Vietoris–Rips Persistence Barcodes," *JACT* (2021)
