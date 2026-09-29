# A Plane-Wave Integral Involving Associated Legendre Functions: Reopening a Cold Case with Modern AI

A previously unresolved problem from my PhD is revisited and explored with Large Language Model (LLM) assistance to (i) revisit the original mathematics and (ii) develop and benchmark a robust numerical scheme. The work is motivated by the elegance found in the original approach, and the authors frustration (at the time) at being unable to utilize it as hoped for.

---

## Repository Structure

- `data/raw/`: Spreadsheet data (`.xlsx` and `.csv`) containing high-precision benchmark evaluations computed via `mpmath` and numerical quadrature comparisons across angular momentum configurations.
- `data/figures/`: CSV data and Python scripts used to generate image files.
- `src/mavromatis_cauchy/`: Core package containing series formulations using `mpmath` and incorporating analytical truncation bounds $k_{\max}(a) = \lceil a + 5a^{1/3} + 15 \rceil$.

---

## Installation

### Prerequisites
- Python 3.10 or later
- Recommended: Virtual environment (`venv` or `conda`)

### Setup with `pip`
```bash
git clone [https://github.com/](https://github.com/)<your-username>/mavromatis-cold-case.git
cd mavromatis-cold-case
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
