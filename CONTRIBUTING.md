# Contributing to Celiac Trial Simulator

Celiac Trial Simulator is an open-source biostatistical framework analyzing whether published Phase 2 celiac disease clinical trials had adequate statistical power to detect efficacious drug candidates.

## Development Setup

The project uses Python 3.12+ and `uv`:

```bash
# Clone
git clone https://github.com/Ahmet-Dedeler/celiac-trial-sim.git
cd celiac-trial-sim

# Install dependencies and sync virtual environment
uv sync

# Run tests
uv run pytest
```

## How to Contribute

- Browse [open issues](https://github.com/Ahmet-Dedeler/celiac-trial-sim/issues) for starter tasks.
- Contributions welcomed for adding newly published celiac challenge trial data, Monte Carlo simulation optimizations, and interactive visualizations.
