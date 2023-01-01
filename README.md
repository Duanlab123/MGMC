The geometric high-order regularization methods such as mean curvature and Gaussian curvature, have been intensively studied during the last decades due to their abilities in preserving geometric properties including image edges, corners, and image contrast. However, the dilemma between restoration quality and computational efficiency is an essential roadblock for high-order methods. In this paper, we propose fast multi-grid algorithms for minimizing both mean curvature and Gaussian curvature energy functionals without sacrificing accuracy for efficiency. Unlike the existing approaches based on operator splitting and the Augmented Lagrangian method (ALM), no artificial parameters are introduced in our formulation, which guarantees the robustness of the proposed algorithm. Meanwhile, we adopt the domain decomposition method to promote parallel computing and use the fine-to-coarse structure to accelerate convergence. Numerical experiments are presented on image denoising, CT, and MRI reconstruction problems to demonstrate the superiority of our method in preserving geometric structures and fine details. The proposed method is also shown effective in dealing with large-scale image processing problems by recovering an image of size  within $40$s, while the ALM method requires around $200$s.
please cite this paper: @article{zhang2022fast,

title={Fast Multi-grid Methods for Minimizing Curvature Energy},

author={Zhang, Zhenwei and Chen, Ke and Duan, Yuping},

journal={arXiv preprint arXiv:2204.07921},

year={2022}

}
