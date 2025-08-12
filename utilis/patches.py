# utils/patches.py
import torch
import collections
import types as _types
import sys as _sys

def apply_patches():
    """Applies all necessary compatibility patches."""
    # FAKE torch.float8 stub (for PyTorch <2.1)
    if not hasattr(torch, "float8_e4m3fn"):
        torch.float8_e4m3fn = torch.float16
    if not hasattr(torch, "float8_e5m2"):
        torch.float8_e5m2 = torch.float16

    # torch.compiler stub for flex_attention (PyTorch 2.0)
    if not hasattr(torch, "compiler"):
        compiler_stub = _types.SimpleNamespace()
        def disable(*args, **kwargs):
            def decorator(fn): return fn
        compiler_stub.disable = disable
        torch.compiler = compiler_stub
    
    # Patch torch._six for legacy compatibility
    try:
        import torch._six as torch_six
    except ImportError:
        torch_six = _types.ModuleType("torch._six")
        _sys.modules["torch._six"] = torch_six
    additions = [
        ("inf", float("inf")),
        ("string_classes", (str,)),
        ("integer_classes", (int,)),
        ("class_types", (type,)),
        ("container_abcs", collections.abc),
        ("moves", _types.SimpleNamespace(collections_abc=collections.abc)),
    ]
    for attr, val in additions:
        if not hasattr(torch_six, attr):
            setattr(torch_six, attr, val)
