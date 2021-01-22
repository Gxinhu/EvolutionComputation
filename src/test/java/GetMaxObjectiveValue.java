import cn.hutool.core.util.StrUtil;
import cn.hutool.json.JSONArray;
import cn.hutool.json.JSONObject;
import cn.hutool.json.JSONUtil;
import jmetal.qualityIndicator.util.MetricsUtil;
import org.apache.commons.math3.linear.MatrixUtils;
import org.junit.jupiter.api.Test;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.toList;

public class GetMaxObjectiveValue {
    public static void main(String[] args) {
        String pfPath = "./10000PF/";
        List<Integer> objectiveList = Arrays.asList(5, 8, 10, 15);
        IntStream.rangeClosed(1, 7)
                .parallel()
                .forEach(
                        i -> objectiveList.stream()
                                .parallel()
                                .map(objective ->
                                        StrUtil.format("{}MaF/{}d/MaF{}.pf", pfPath, objective, i))
                                .map(MetricsUtil::readFront)
                                .filter(Objects::nonNull)
                                .map(MatrixUtils::createRealMatrix)
                                .map(realMatrix -> IntStream.range(0, realMatrix.getColumnDimension()).
                                        parallel()
                                        .mapToObj(j -> realMatrix.getColumnVector(j).getMaxValue())
                                        .collect(toList()))
                                .map(doubles -> StrUtil.format("MaF{}-{}d:{}", i, doubles.size(), doubles))
                                .sorted()
//                                .forEach(System.out::println)
                );

    }

    @Test
    public void testJSON() {

        int numberOfObjective = 5;

        String problemNames = "MaF1";
        if (problemNames.contains("MaF")) {
            File path = new File("prop/maxObjectiveValue.json");
            JSONObject maxPoint = JSONUtil.readJSONObject(path, StandardCharsets.UTF_8);
            String fullName = StrUtil.format("{}-{}d", problemNames, numberOfObjective);
            JSONArray points = maxPoint.getJSONArray(fullName);
        }
    }
}
