# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['ISM_App.py'],
             pathex=['C:\\Users\\Mo\\Dropbox\\notes\\Bioinformatics\\Motif_Search_GUI'],
             binaries=[],
             datas=[('C:\\Users\\Mo\\Dropbox\\notes\\Bioinformatics\\Motif_Search_GUI\\MM_icon.ico', '.')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='ISM_App',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False , icon='apta_index.ico')
