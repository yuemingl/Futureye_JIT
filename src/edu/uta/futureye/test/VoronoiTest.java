package edu.uta.futureye.test;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.GridLayout;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.UIManager;
import javax.swing.event.MouseInputAdapter;


public class VoronoiTest extends MouseInputAdapter {
	
	
	JPanel pane = null;
	public Component createComponents() {
		final JLabel label = new JLabel("局部坐标 (r,s) => 物理坐标 (x,y)");
		final JButton btn = new JButton("Start");
		pane = new JPanel();
		pane.setBorder(BorderFactory.createEmptyBorder(50, // top
				50, // left
				600, // bottom
				600) // right
				);
		pane.setLayout(new GridLayout(0, 1)); // 单列多行
		pane.add(label);
		pane.add(btn);
		
		//pane.addMouseMotionListener(this);
		pane.addMouseListener(this);	
		return pane;
	}
	public void mouseMoved(MouseEvent e) {
		int x = e.getX();
		int y = e.getY();
		pane.getGraphics().drawOval(x, y, 5, 5);
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub
		int x = e.getX();
		int y = e.getY();
		pane.getGraphics().drawOval(x, y, 5, 5);
		
	}

	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		try {
			UIManager.setLookAndFeel(UIManager
					.getCrossPlatformLookAndFeelClassName());
			// 设置窗口风格
		} catch (Exception e) {
		}

		// 创建顶层容器并添加内容.
		JFrame frame = new JFrame("Voronoi Test");
		VoronoiTest app = new VoronoiTest();
		Component contents = app.createComponents();
		frame.getContentPane().add(contents, BorderLayout.CENTER);

		// 窗口设置结束，开始显示
		frame.addWindowListener(new WindowAdapter() {
			// 匿名类用于注册监听器
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		frame.pack();
		frame.setVisible(true);
	}
}
