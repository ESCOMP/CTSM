.. _running-with-ftorch-section:

###############################################
Running with Machine Learning Models via FTorch
###############################################

FTorch (Fortran interface to PyTorch) allows CTSM Fortran routines to call PyTorch models at runtime for machine learning parameterizations, emulators, and data-driven subgrid processes.

Prerequisites
=============

1. **Python Environment (for model export)**:
   PyTorch is only required if you are creating or exporting models into TorchScript format (CTSM runtime simulations link directly against C++ LibTorch and do not require PyTorch in your Python environment). To export models, install PyTorch into your active Python environment:

   * In ``ctsm_pylib`` or another conda environment:

     .. code-block:: bash

        conda activate ctsm_pylib
        conda install pytorch-cpu -c pytorch
        # or via pip:
        pip install torch --index-url https://download.pytorch.org/whl/cpu

   * On NCAR HPC systems (e.g., Derecho), you can alternatively activate NCAR's pre-built Python environment with PyTorch pre-installed:

     .. code-block:: bash

        module load conda
        conda activate npl

2. **Checking out FTorch**:
   Ensure the ``libraries/FTorch`` submodule is populated via ``git-fleximod``:

   .. code-block:: bash

      ./bin/git-fleximod update

Workflow: Exporting PyTorch Models to TorchScript (.pt)
=======================================================

FTorch requires models to be exported into TorchScript format (``.pt``). You can export custom PyTorch ``nn.Module`` models using ``torch.jit.script`` / ``torch.jit.trace`` (or ``ctsm.machine_learning.export_torchscript.export_model_to_torchscript``), or generate a built-in demonstration model (``constant`` or ``identity``) using the provided CTSM export tool:

.. code-block:: bash

   ./tools/ftorch/export_torchscript -o my_model.pt

At runtime, ensure the exported ``.pt`` model file is staged into your case's run directory (``$RUNDIR``).

Building CTSM with FTorch Support
=================================

To build CTSM with FTorch enabled, set ``USE_FTORCH=TRUE`` in your case configuration before building:

.. code-block:: bash

   ./xmlchange USE_FTORCH=TRUE
   ./case.setup
   ./case.build

Interfacing Fortran Code with FTorch
====================================

To load and execute a PyTorch model in CTSM Fortran code, import from ``FTorch_cesm_interface`` (see ``use_ftorch_demo`` for a working example):

.. code-block:: fortran

   use FTorch_cesm_interface, only : torch_model, torch_tensor, torch_model_load, &
                                     torch_tensor_from_array, torch_model_forward, &
                                     torch_kCPU, torch_delete

   type(torch_model) :: model
   type(torch_tensor) :: in_tensor(1), out_tensor(1)
   real(r8), target :: input_data(N), output_data(N)
   integer :: in_layout(1), out_layout(1)

   in_layout = [1]
   out_layout = [1]

   ! 1. Load TorchScript model file (from $RUNDIR)
   ! (do this once, in an initialize step.)
   call torch_model_load(model, "my_model.pt", torch_kCPU)

   ! 2. Bind Fortran arrays to Torch tensors (shape is inferred from array; layout orders dimensions)
   ! (do this and #3 once per inference - probably once per timestep)
   call torch_tensor_from_array(in_tensor(1), input_data, in_layout, torch_kCPU)
   call torch_tensor_from_array(out_tensor(1), output_data, out_layout, torch_kCPU)

   ! 3. Perform forward inference pass
   call torch_model_forward(model, in_tensor, out_tensor)

   ! 4. Clean up model and tensor handles
   ! (do this once, in a teardown step.)
   call torch_delete(in_tensor(1))
   call torch_delete(out_tensor(1))
   call torch_delete(model)
