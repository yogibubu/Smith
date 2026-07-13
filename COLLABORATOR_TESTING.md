# Private collaborator test for standalone SMITH

The purpose of this first external test is to find installation, portability,
command-line, and documentation problems before deciding whether the manuscript
needs additional scientific examples.  The full ORACLE application is still in
final testing and is not part of this test.

## 1. Give the collaborator access

Keep `yogibubu/Smith` private.  In GitHub open:

`Settings` -> `Collaborators` -> `Add people`

Enter the collaborator's GitHub username.  They must accept the invitation
before cloning or installing the private repository.

## 2. Clone and install

The collaborator should install GitHub CLI, then run:

```bash
gh auth login
gh auth setup-git
git clone https://github.com/yogibubu/Smith.git
cd Smith
git switch agent/oracle-boundary-standalone-smith
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install ./standalone
```

On Windows PowerShell, activation is:

```powershell
.venv\Scripts\Activate.ps1
```

After this test branch is merged, the `git switch` line is no longer needed.

## 3. Run the packaged examples

```bash
smith-sonic --version
smith-sonic example water water.xyzin
smith-sonic inspect water.xyzin
smith-sonic example norbornane norbornane.xyzin
smith-sonic inspect norbornane.xyzin
```

Expected summaries are:

- water: 3 GICs, rank 3;
- norbornane: 51 GICs, rank 51;
- both outputs: `PERCEPTION_PROFILE REDUCED_ORACLE`.

The source inputs are also visible under `standalone/examples/`.

## 4. Information to return

Please report:

- operating system, processor architecture, and Python version;
- whether authentication, cloning, and installation completed without manual
  intervention;
- the complete terminal output if a command fails;
- whether the distinction between reduced ORACLE perception and SMITH/SONIC
  coordinate construction is understandable;
- whether the commands and generated files are self-explanatory;
- the summaries printed for water and norbornane.

Do not broaden this first test to the unreleased full ORACLE application.  A
non-covalent complex and an eta3 transition-metal complex are candidate SMITH
examples to evaluate only after this basic external test has identified any
packaging or interface problems.

