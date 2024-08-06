# &infin;tools 🧰
👷‍♂️ Work in progress python package that provide various AIRETIS tools for tasks like setting up simulations and analyzing data. 📈

🤝 Contributions of tools and modifications are welcome – feel free to submit your enhancements! 🚀
## Installation 

1. Clone the repository
2. (Optional) Create a virtual environment
3. Install with:
   ```bash   
   python -m pip install -e .
   ```
## Usage
```bash
inft -h
```

## Development practice

Current CLI capabilities is powered by the [Typer](https://typer.tiangolo.com/) python library.

Any function `function_name` added to the folders [`inftools/exercises`](`inftools/exercises`) and [`inftools/tistools`](`inftools/tistools`) will automatically be callable via the terminal as

```bash
inft function_name -args1 var1 -args2 var2
```

The function should have the following structure:

```python3
from typing import Annotated

def function_name(
    args1: Annotated[str, typer.Option("-args1", help="var1")],
    args2: Annotated[str, typer.Option("-args2", help="var2")],
    ):
    """
    Function description
    """
```

Individual functions can also be added to the CLI library by including them in [`inftools/bin.py`], as done for the wham function in the [`inftools/analysis`](`inftools/analysis`) folder:

```python3
from inftools.analysis.wham import wham

MAPPER["wham"] = wham
```
