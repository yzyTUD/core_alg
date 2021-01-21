package algo_paralleltest_helloworld;

public class Counter {
	public static int val=0;
	static synchronized void incrementByOne() {
		val+=1;
	}
}


package algo_paralleltest_helloworld;

public class FirstThreadTest extends Thread {
	
    public void run() {
        for (int i = 0; i < 100; i++) {
            //System.out.println(this.getName() + " " + i);
        	Counter.incrementByOne();
        }
    }

    public static void main(String[] args) {
    	FirstThreadTest ft = new FirstThreadTest();
    	ft.start();
        for (int i = 0; i < 100; i++) {
            //System.out.println(Thread.currentThread().getName() + " " + i);
        	Counter.incrementByOne();
        }
        
        try {
			ft.join();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        System.out.println(Counter.val);
    }
}