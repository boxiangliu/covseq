from setuptools import setup

setup(
    name='coviz',
    version='0.1',
    py_modules=['coviz'],
    install_requires=[
        'Click',
    ],
    entry_points='''
        [console_scripts]
        annotate=annotation:annotate
        predict=orf_prediction:predict
    ''',
)