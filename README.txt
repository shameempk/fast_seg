Dependencies:

1. System wide dependencies
	1. Python 3.5+
	2. opencv (along with python bindings, comes together by default)
2. Project dependencies
	(Recommends using python virtualenv: http://docs.python-guide.org/en/latest/dev/virtualenvs )
	1. All dependencies are enlisted in requirements.txt
		Install them using : pip install -r requirements.txt


Executing the code:

1. Run the main file using python3: python3 fast_seg.py
	* Will provide a minimal GUI to mark the seed pixels. While marking, switching between "background" and "object" pixels are done using keys 'b' and 'o' respectively. By default GUI initializes in object mode. Object is marked with "red" markings and Background with "blue".
2. Press ESC after marking the seeds.
3. Output window will provide the results.
4. Output image will be written in running folder, named "out.png"


For any other inquiries file an issue at https://github.com/shameempk/fast_seg .