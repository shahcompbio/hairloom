import subprocess
import importlib.resources as pkg_resources


def run_cli_command(command):
    """Run a CLI command and return the output."""
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return result.stdout.strip(), result.stderr.strip(), result.returncode


def load_expected_output(file_path):
    """Load the expected output from a file."""
    with open(file_path, "r") as f:
        return f.read().strip()


def test_cli_extract():
    """Test the 'breakpoints' CLI command."""
    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test_reads.bam'))
    expected_file = str(pkg_resources.files("hairloom").joinpath('data/cli/extract.txt'))
    command = ["hairloom", "extract", bam_path, "chr1", "50", "150"]

    stdout, stderr, returncode = run_cli_command(command)
    expected_output = load_expected_output(expected_file)

    assert returncode == 0, f"CLI failed with stderr: {stderr}"
    assert stdout == expected_output, f"Output does not match expected output:\n{stdout}"


def test_cli_breakpoints():
    """Test the 'breakpoints' CLI command."""
    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test_reads.bam'))
    expected_file = str(pkg_resources.files("hairloom").joinpath('data/cli/breakpoints.txt'))
    command = ["hairloom", "breakpoints", bam_path, "chr1", "50", "150"]

    stdout, stderr, returncode = run_cli_command(command)
    expected_output = load_expected_output(expected_file)

    assert returncode == 0, f"CLI failed with stderr: {stderr}"
    assert stdout == expected_output, f"Output does not match expected output:\n{stdout}"


def test_cli_segments():
    """Test the 'segments' CLI command."""
    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test_reads.bam'))
    expected_file = str(pkg_resources.files("hairloom").joinpath('data/cli/segments.txt'))
    command = ["hairloom", "segments", bam_path, "chr1", "50", "150"]

    stdout, stderr, returncode = run_cli_command(command)
    expected_output = load_expected_output(expected_file)

    assert returncode == 0, f"CLI failed with stderr: {stderr}"
    assert stdout == expected_output, f"Output does not match expected output:\n{stdout}"


def test_cli_svs():
    """Test the 'svs' CLI command."""
    bam_path = str(pkg_resources.files("hairloom").joinpath('data/test_reads.bam'))
    expected_file = str(pkg_resources.files("hairloom").joinpath('data/cli/svs.txt'))
    command = ["hairloom", "svs", bam_path, "chr1", "50", "150"]

    stdout, stderr, returncode = run_cli_command(command)
    expected_output = load_expected_output(expected_file)

    assert returncode == 0, f"CLI failed with stderr: {stderr}"
    assert stdout == expected_output, f"Output does not match expected output:\n{stdout}"