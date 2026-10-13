import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load_module(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class ContainerRuntimeTests(unittest.TestCase):
    def test_manifest_consistency(self):
        m = json.loads((ROOT / 'bin/containers/images.json').read_text())
        for tool, details in m.items():
            self.assertTrue(details['docker'])
            self.assertTrue(details['apptainer'])
            if details['dockerfile']:
                self.assertTrue((ROOT / details['dockerfile']).exists(), tool)
            else:
                self.assertEqual(tool, 'busco')
            if details['definition']:
                self.assertTrue((ROOT / details['definition']).exists(), tool)
            else:
                self.assertEqual(tool, 'busco')

    def test_docker_command_mounts_and_ownership(self):
        mod = load_module('docker_exec', 'bin/scripts/docker_exec.py')
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'repo'; root.mkdir()
            outside = Path(temp) / 'reads'; outside.mkdir()
            command = mod.docker_command('image:tag', root, [str(outside), str(root)], 'echo Hello')
            self.assertEqual(command[:3], ['docker', 'run', '--rm'])
            self.assertIn('--user', command)
            self.assertIn(f'{os.getuid()}:{os.getgid()}', command)
            self.assertIn(f'type=bind,src={root},dst={root}', command)
            self.assertIn(f'type=bind,src={outside},dst={outside},readonly', command)
            self.assertEqual(command[-1], 'echo Hello')
            self.assertIn('--entrypoint', command)
            self.assertEqual(command[command.index('--entrypoint')+1], '/bin/bash')

    def test_docker_command_rejects_missing_mount(self):
        mod = load_module('docker_exec_other', 'bin/scripts/docker_exec.py')
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileNotFoundError):
                mod.docker_command('x:y', Path(temp), ['/missing/does/not/exist'], 'true')

    def test_launcher_dry_run_docker_and_apptainer(self):
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            db = work / 'nr.dmnd'; db.write_text('fake')
            fq = work / 'pr_197_001.hifireads.fastq.gz'; fq.write_text('fake')
            bin_dir = work / 'bin'; bin_dir.mkdir()
            captured = work / 'args.json'
            snk = bin_dir / 'snakemake'
            snk.write_text('#!/usr/bin/env python3\nimport json,sys,os\njson.dump(sys.argv[1:],open(os.environ["SNK_CAPTURE"],"w"))\n')
            snk.chmod(0o755)
            env = os.environ.copy()
            env['PATH'] = str(bin_dir) + os.pathsep + env['PATH']
            env['SNK_CAPTURE'] = str(captured)
            for runtime in ('docker', 'apptainer', 'none'):
                cmd = ['bash', 'run_pipeline.sh', '--reads', str(fq), '--db', str(db),
                       '--container-runtime', runtime, '--dry-run']
                result = subprocess.run(cmd, cwd=ROOT, env=env, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                args = json.loads(captured.read_text())
                self.assertIn('-n', args)
                self.assertIn('container_runtime=' + runtime, args)
                self.assertEqual('--use-apptainer' in args, runtime == 'apptainer')
                self.assertIn('pr_197_001', result.stdout)
                self.assertNotIn('setup_containers.py', result.stdout)

    def test_launcher_standalone_dry_run_does_not_require_sif(self):
        with tempfile.TemporaryDirectory() as td:
            work = Path(td)
            db = work / 'nr.dmnd'; db.touch()
            fq = work / 's.fastq.gz'; fq.write_text('fake')
            ref = work / 'ref.faa'; ref.touch()
            bindir = work / 'bin'; bindir.mkdir()
            mock = bindir / 'snakemake'
            mock.write_text('#!/bin/sh\nexit 0\n'); mock.chmod(0o755)
            env = os.environ.copy(); env['PATH'] = str(bindir) + os.pathsep + env['PATH']
            result = subprocess.run(['bash','run_pipeline.sh','--reads',str(fq),'--db',str(db),
                                     '--gene-prediction','standalone','--metaeuk-ref',str(ref),
                                     '--container-runtime','docker','--dry-run'],
                                    cwd=ROOT,env=env,text=True,capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_setup_uses_cached_images(self):
        mod = load_module('setup_containers','bin/scripts/setup_containers.py')
        manifest = json.loads((ROOT / 'bin/containers/images.json').read_text())
        from unittest.mock import patch
        with patch.object(mod.shutil, 'which', return_value='/usr/bin/docker'), \
             patch.object(mod.subprocess, 'run') as run:
            run.return_value.returncode = 0
            mod.provision(ROOT, 'docker', manifest, ['longqc'], False)
            self.assertEqual(run.call_count, 2)  # docker info, docker image inspect

    def test_snakefile_uses_native_or_docker_commands(self):
        s = (ROOT / 'bin/Snakefile').read_text()
        for tool in ('longqc','flye','coverm','busco','telfinder','diamond','report','metaeuk','prodigal','tiara'):
            self.assertIn(f'container_image("{tool}")', s)
            self.assertIn(f'tool_command("{tool}", ', s)
        self.assertNotIn('REPORT_SIF=', s)


if __name__ == '__main__': unittest.main()
