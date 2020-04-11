<h1>Open Snakes for strip like object segmentation</h1>
Snakes for segmentation of strip-like object (in this work, we take intima media borders in ultrasound sequences as an example) with different morphologies and dynamics.
<h2>Brief Introduction: </h2>
This repository is the MATLAB implementation of paper "<a href="https://scholar.google.com.hk/citations?hl=zh-CN&user=U68XZOgAAAAJ&view_op=list_works&sortby=pubdate" target="_blank"><b> Robust Segmentation of Intima-Media Borders with Different Morphologies and Dynamics During the Cardiac Cycle. </b></a>" (published in IEEE Journal of Biomedical and Health Informatics (JBHI) 2018). <br>
This code provides a robust, accurate and efficient method based on snakes (active contour models) of segmenting Intima-Media borders from ultrasound sequences.<br>

<h2>Code Usage: </h2>
<h3>Inputs:  </h3>
<ul>
    <li>an ultrasound sequence
    <li>a manual segmentation of the first frame
</ul>
<h3>Outputs:  </h3>
<ul>
    <li>the segmentation of the whole sequence
</ul>
<h3>Usage:  </h3>
Simply download the files, make sure the path to the sequence is correct, and run the "main_snake.m" file in folder "sequence_segmentation_using_manual_initial". <br>
<h3>Other Notes:  </h3>
<ul>
<li>Almost all comments in the codes are written in Chinese now. If I have time in the future, the comments might be translated into English. However, to be frank, the translation may be suspended indefinitely since I get more and more work to do. Fortunately, you can find someone who knows Chinese almost everywhere around the world. :D<br>
<li>To make the codes applicable to sequences of more morphologies and dynamics, I have edited the codes after the paper is published. Thus, you may find some details in the codes that are different from the paper. However, the idea is the same, i.e, the external energy matrix of the snake is modified based on the grayscale and derivatives.<br>
<li>A MATLAB program for automatically obtaining the segmentation of the first frame will be made on-line soon, which can replace the manual segmented IM borders for the first frame.
<li>All questions and debuggings are welcome! If you find the code useful, please kindly cite my paper (BibTeX citations are as follows).<br>
	@article{zhao2017robust,
	  title={Robust segmentation of intima--media borders with different morphologies and dynamics during the cardiac cycle},
	  author={Zhao, Shen and Gao, Zhifan and Zhang, Heye and Xie, Yaoqin and Luo, Jianwen and Ghista, Dhanjoo and Wei, Zhanghong and Bi, Xiaojun and Xiong, Huahua and Xu, Chenchu and others},
	  journal={IEEE journal of biomedical and health informatics},
	  volume={22},
	  number={5},
	  pages={1571--1582},
	  year={2017},
	  publisher={IEEE}
	}
</ul>