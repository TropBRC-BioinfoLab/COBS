from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional, Literal

ParamType = Literal["int", "float", "bool", "select"]

@dataclass
class ParamSpec:
    key: str
    label: str
    ptype: ParamType
    default: Any
    help: str = ""

    # slider untuk int/float
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    step: Optional[float] = None

    # untuk select
    options: Optional[List[Any]] = None

    def render(self, st) -> Any:
        if self.ptype == "bool":
            return st.checkbox(self.label, value=bool(self.default), help=self.help)

        if self.ptype == "select":
            opts = self.options or []
            idx = 0
            if self.default in opts:
                idx = opts.index(self.default)
            return st.selectbox(self.label, opts, index=idx, help=self.help)

        if self.ptype == "int":
            return st.slider(
                self.label,
                int(self.min_value if self.min_value is not None else 0),
                int(self.max_value if self.max_value is not None else 100),
                int(self.default),
                int(self.step if self.step is not None else 1),
                help=self.help,
            )

        if self.ptype == "float":
            return st.slider(
                self.label,
                float(self.min_value if self.min_value is not None else 0.0),
                float(self.max_value if self.max_value is not None else 1.0),
                float(self.default),
                float(self.step if self.step is not None else 0.05),
                help=self.help,
            )

        raise ValueError(f"Unknown ptype: {self.ptype}")

@dataclass
class AlgoSpec:
    key: str
    name: str
    description: str
    params: List[ParamSpec]
    run: Callable  # run(G, params_dict) -> NodeClustering
